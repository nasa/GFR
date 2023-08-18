module metrics_mod
  !
!#define TIME_METRIC_PROCEDURES
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
  type, public :: face_nrm_t
   !real(wp) :: area
    real(wp), allocatable :: v(:,:)
  contains
    procedure, public :: area => compute_face_area
  end type face_nrm_t
  !
  ! geofa : normal of the face at each of the solution points on the face
  !
  type(face_nrm_t), public, save, allocatable :: geofa(:)
  !
  ! facenrm : normal of the face at each of the flux points
  !
  real(wp), public, save, allocatable :: facenrm(:,:)
  real(wp), public, save, allocatable :: facejac(:,:)
  !
  !$ACC DECLARE CREATE(facenrm,facejac)
  !
  ! facenrm_qdtpts : normal of the face at each of the face quadrature points
  !
  real(wp), public, save, allocatable :: facenrm_qdtpts(:,:)
  !
  !$ACC DECLARE CREATE(facenrm_qdtpts)
  !
  ! face_rrf_vel : array containing the precomputed value of
  !                      crossproduct(omega,xyz)
  !                at each flux point where
  !                   omega = the rotation rate vector of the
  !                           rotating reference frame
  !                   xyz = Cartesian coordinates at the flux point
  ! face_rrf_vel_qdtpts : same as face_rrf_vel but for face quadrature points
  ! NOTE: These two arrays have been created because the angular velocity
  !       is needed at each face/flux point (FP) for a rotating reference frame
  !       in order to compute the interface fluxes. The coordinates of each FP
  !       is needed to compute this angular velocity, and these coordinates
  !       arent usually saved to an array. Instead of saving these coordinates
  !       to the facexyz array and then recomputing the angular velocity at each
  !       FP each time it is needed, we are instead just precomputing the
  !       angular velocity and storing that.
  !
  real(wp), public, save, allocatable :: face_rrf_vel(:,:)
  real(wp), public, save, allocatable :: face_rrf_vel_qdtpts(:,:)
  !
  !$ACC DECLARE CREATE(face_rrf_vel,face_rrf_vel_qdtpts)
  !
  ! wtdjac : jacobian * quadrature weight at each solution point
  !
  real(wp), public, save, allocatable :: wtdjac(:)
  !
  ! jac : jacobian at each solution point
  !
  real(wp), public, save, allocatable :: jac(:)
  !
  !$ACC DECLARE CREATE(wtdjac,jac)
  !
  ! time_metrics : delta_x, delta_y, and delta_z components used to compute
  !                the point-wise time-step using the structured grid CFL
  !                condition.
  !
  !       dt = 1 / [ abs(u)/dx + abs(v)/dy + abs(w)/dz + &
  !                  a*sqrt(1/(dx**2) + 1/(dy**2) + 1/(dz**2)) ]
  !       where:
  !         u = x-velocity
  !         v = y-velocity
  !         w = z-velocity
  !         a = speed of sound [a=sqrt(gam*pressure/density) ]
  !         time_metrics(1,:) = 1/dx
  !         time_metrics(2,:) = 1/dy
  !         time_metrics(3,:) = 1/dz
  !         time_metrics(4,:) = sqrt(1/(dx**2) + 1/(dy**2) + 1/(dz**2))
  !
  real(wp), public, save, allocatable :: time_metrics(:,:)
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
  real(wp), public, save, allocatable :: metrics_array(:,:,:)
  !$ACC DECLARE CREATE(metrics_array)
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
  use geovar, only : nr,ncell,nnode
  use geovar, only : nfbnd,nface,n_solpts
  use geovar, only : using_quadrature
  use ovar,   only : local_volume,total_volume
  use ovar,   only : Timestep_Type
  !
  !.. Local Scalars ..
  integer  :: ierr
  integer  :: min_jsp,min_jcl,min_jac_cpu,min_vol_cpu
  integer  :: max_jsp,max_jcl,max_jac_cpu,max_vol_cpu
  integer  :: min_ncell_cpu,max_ncell_cpu
  integer  :: min_nnode_cpu,max_nnode_cpu
  integer  :: min_nfbnd_cpu,max_nfbnd_cpu
  integer  :: min_nface_cpu,max_nface_cpu
  integer  :: min_nspts_cpu,max_nspts_cpu
  integer  :: min_ncv,max_ncv,min_ncv_cpu,max_ncv_cpu
  real(wp) :: ave_ncell,ave_nnode,ave_vol,ave_cell_vol
  real(wp) :: ave_nfbnd,ave_nface,ave_nspts
  real(wp) :: jac_min,jac_max,vol_min,vol_max
  !
  !.. Local Arrays ..
  integer  :: part_info(1:5)
  integer  :: all_part_info(1:5,1:ncpu)
  real(wp) :: jac_info(1:11)
  real(wp) :: all_jac_info(1:11,1:ncpu)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_metrics"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  jac_info(:) = zero
  all_jac_info(:,:) = zero
  !
  ! Allocate the required metrics related arrays
  !
  allocate ( wtdjac(1:n_solpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"wtdjac",1,__LINE__,__FILE__,ierr,error_message)
  !$ACC UPDATE DEVICE(wtdjac)
  !
  allocate ( jac(1:n_solpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"jac",1,__LINE__,__FILE__,ierr,error_message)
  !$ACC UPDATE DEVICE(jac)
  !
  allocate ( metrics_array(1:n_solpts,1:nr,1:nr) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"metrics_array",1,__LINE__,__FILE__,ierr,error_message)
  !$ACC UPDATE DEVICE(metrics_array)
  !
  !
  ! Loop over the cells to get the metrics at each solution point
  !
  if (nr == 3) then
    call compute_cell_metrics_3d(metrics_array,jac,wtdjac)
  else if (nr == 2) then
    call compute_cell_metrics_2d(metrics_array,jac,wtdjac)
  else if (nr == 1) then
    call compute_cell_metrics_1d(metrics_array,jac,wtdjac)
  else
    write (error_message,100) nr
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !$ACC UPDATE DEVICE(jac,wtdjac,metrics_array)
  !
  ! Compute volume and jacobian statistics
  !
  call compute_volume_statistics(min_ncv,max_ncv,vol_min,vol_max, &
                                 min_jcl,max_jcl,min_jsp,max_jsp, &
                                 jac_min,jac_max,jac,wtdjac)
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
  ! Get the total volume of the grid by summing all the local partition volumes
  !
  total_volume = sum( all_jac_info(7,:) )
  !
  ! Have the global root process its value of the total volume to all other
  ! processes just to make sure all processes have the same value
  !
  if (ncpu > 1) then
    call mpi_bcast(total_volume,1_int_mpi,mpi_flttyp, &
                   glb_root,MPI_COMM_WORLD,mpierr)
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
  if (mypnum == glb_root) then
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
  ! Get the face normal at the solution points on each face
  !
  call compute_face_normals
  !
  ! Get the face normal at the flux quadrature points on each face
  !
  if (using_quadrature) then
    call compute_quadrature_face_normals
  end if
  !
  if (abs(Timestep_Type) == Point_Timestep) then
    call create_time_metrics
  end if
  !
  call mpi_barrier(MPI_COMM_WORLD,mpierr)
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
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
  100  format ("ERROR in subroutine get_metrics! The value of nr ", &
               "should be in the range of [1-3], but nr = ",i0,"!")
  !
end subroutine get_metrics
!
!###############################################################################
!
subroutine old_get_metrics()
  !
  !.. Use Statements ..
  use geovar, only : nr,n_totpts,n_totcel
  use geovar, only : ncell,cell,xyz
  use geovar, only : nnode,nfbnd,nface,n_solpts
  use geovar, only : using_quadrature
  !
  use element_type_mod, only : std_elem
  !
  use ovar, only : local_volume,total_volume
  use ovar, only : VortexGrid_random_perturb
  use ovar, only : Timestep_Type
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
  character(len=*), parameter :: pname = "old_get_metrics"
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
#ifdef TIME_METRIC_PROCEDURES
  call debug_timer(reset_accumulated_time,"output metrics",acc=9)
  call debug_timer(reset_accumulated_time,"setup and transposing xyz",acc=2)
  call debug_timer(reset_accumulated_time,"Filling individual metrics_dt(nc)", &
                   acc=3)
  call debug_timer(reset_accumulated_time, &
                   "Migrating metrics_dt to jac, wtdjac, and metrics_array", &
                   acc=4)
  call debug_timer(reset_accumulated_time,"computing cell volumes",acc=5)
  call debug_timer(reset_accumulated_time,"computing volume statistics",acc=6)
#endif
  !
#ifdef TIME_METRIC_PROCEDURES
  call debug_timer(start_timer,"Metric initialization operations")
#endif
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
#ifdef TIME_METRIC_PROCEDURES
    call debug_timer(start_accumulation,"output metrics",acc=9)
#endif
    do i = ibeg,iend
      open(newunit=iom(i),file=fname(i),status="replace",action="write", &
                          form="formatted",position="rewind", &
                          iostat=ierr,iomsg=error_message)
      call io_error(pname,fname(i),1,__LINE__,__FILE__,ierr,error_message)
    end do
#ifdef TIME_METRIC_PROCEDURES
    call debug_timer(stop_accumulation,"output metrics",acc=9)
#endif
  end if
  !
  ! Allocate the required metrics related arrays
  !
  allocate ( metrics_dt(1:ncell) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"metrics_dt",1,__LINE__,__FILE__,ierr,error_message)
  !
  allocate ( wtdjac(1:n_solpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"wtdjac",1,__LINE__,__FILE__,ierr,error_message)
  !$ACC UPDATE DEVICE(wtdjac)
  !
  allocate ( jac(1:n_solpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"jac",1,__LINE__,__FILE__,ierr,error_message)
  !$ACC UPDATE DEVICE(jac)
  !
  allocate ( metrics_array(1:n_solpts,1:nr,1:nr) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"metrics_array",1,__LINE__,__FILE__,ierr,error_message)
  !$ACC UPDATE DEVICE(metrics_array)
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
#ifdef TIME_METRIC_PROCEDURES
  call debug_timer(stop_timer,"Metric initialization operations")
#endif
  !
  ! Loop over the cells to get the metrics at each solution point
  !
#ifdef TIME_METRIC_PROCEDURES
  call debug_timer(start_timer,"Main metrics loop")
#endif
  metrics_loop: do nc = 1,ncell
    !
#ifdef TIME_METRIC_PROCEDURES
    call debug_timer(start_accumulation,"setup and transposing xyz",acc=2)
#endif
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    np = n2-n1+1
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Create transpose of the coordinates for this cell to reduce cache misses
    !
    do k = n1,n2
      do l = 1,nr
        txyz(k-n1+1,l) = xyz(l,k)
      end do
    end do
#ifdef TIME_METRIC_PROCEDURES
    call debug_timer(stop_accumulation,"setup and transposing xyz",acc=2)
#endif
    !
    if (output_metrics) then
#ifdef TIME_METRIC_PROCEDURES
      call debug_timer(start_accumulation,"output metrics",acc=9)
#endif
      do i = ibeg,iend
        write (iom(i),21) nc
      end do
#ifdef TIME_METRIC_PROCEDURES
      call debug_timer(stop_accumulation,"output metrics",acc=9)
#endif
    end if
    !
    ! Compute the metrics and jacobian for this cell
    !
#ifdef TIME_METRIC_PROCEDURES
    call debug_timer(start_accumulation, &
                     "Filling individual metrics_dt(nc)",acc=3)
#endif
    call metrics_dt(nc)%create_cell_metrics(txyz(1:np,1:nr), &
                                            this_geom,this_order,nc)
#ifdef TIME_METRIC_PROCEDURES
    call debug_timer(stop_accumulation, &
                     "Filling individual metrics_dt(nc)",acc=3)
#endif
    !
#ifdef TIME_METRIC_PROCEDURES
    call debug_timer(start_accumulation, &
                     "Migrating metrics_dt to jac, wtdjac, and metrics_array", &
                     acc=4)
#endif
    do k = 1,np
      n = n1-1+k
      jac(n) = metrics_dt(nc)%jac(k)
      wtdjac(n) = metrics_dt(nc)%jac(k) &
                * std_elem(this_geom,this_order)%wts(k)
      do l = 1,nr
        do i = 1,nr
          metrics_array(n,i,l) = metrics_dt(nc)%met(k,i,l)
        end do
      end do
    end do
#ifdef TIME_METRIC_PROCEDURES
    call debug_timer(stop_accumulation, &
                     "Migrating metrics_dt to jac, wtdjac, and metrics_array", &
                     acc=4)
#endif
    !
    ! Sum the product of the Jacobian and quadrature weight for each point
    ! in this cell to compute the volume of the current cell
    !
#ifdef TIME_METRIC_PROCEDURES
    call debug_timer(start_accumulation,"computing cell volumes",acc=5)
#endif
    vol = metrics_dt(nc)%vol(this_geom,this_order)
#ifdef TIME_METRIC_PROCEDURES
    call debug_timer(stop_accumulation,"computing cell volumes",acc=5)
#endif
    !
    ! Add the volume of the current cell to the total volume
    ! of the current local partition
    !
#ifdef TIME_METRIC_PROCEDURES
    call debug_timer(start_accumulation,"computing volume statistics",acc=6)
#endif
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
#ifdef TIME_METRIC_PROCEDURES
    call debug_timer(stop_accumulation,"computing volume statistics",acc=6)
#endif
    !
    ! Output the metrics if needed
    !
    if (output_metrics) then
#ifdef TIME_METRIC_PROCEDURES
      call debug_timer(start_accumulation,"output metrics",acc=9)
#endif
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
#ifdef TIME_METRIC_PROCEDURES
      call debug_timer(stop_accumulation,"output metrics",acc=9)
#endif
    end if
    !
  end do metrics_loop
  !$ACC UPDATE DEVICE(jac,wtdjac,metrics_array)
#ifdef TIME_METRIC_PROCEDURES
  call debug_timer(stop_timer,"Main metrics loop")
#endif
  !
#ifdef TIME_METRIC_PROCEDURES
  call debug_timer(start_timer,"Computing and reporting metric statistics")
#endif
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
#ifdef TIME_METRIC_PROCEDURES
  call debug_timer(stop_timer,"Computing and reporting metric statistics")
#endif
  !
#ifdef TIME_METRIC_PROCEDURES
  call debug_timer(report_accumulated_time,"setup and transposing xyz",acc=2)
  call debug_timer(report_accumulated_time, &
                   "Filling individual metrics_dt(nc)", &
                   acc=3)
  call debug_timer(report_accumulated_time, &
                   "Migrating metrics_dt to jac, wtdjac, and metrics_array", &
                   acc=4)
  call debug_timer(report_accumulated_time,"computing cell volumes",acc=5)
  call debug_timer(report_accumulated_time,"computing volume statistics",acc=6)
#endif
  !
  ! Get the face normal at the solution points on each face
  !
  call old_compute_face_normals
  !
  ! Get the face normal at the flux quadrature points on each face
  !
  if (using_quadrature) then
    call compute_quadrature_face_normals
  end if
  !
  if (abs(Timestep_Type) == Point_Timestep) then
    call create_time_metrics
  end if
  !
  ! Reduce the volumes of each partition to the root
  ! processor to determine the volume of the global grid.
  !
#ifdef TIME_METRIC_PROCEDURES
  call debug_timer(start_timer,"Reducing local and total volumes")
#endif
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
#ifdef TIME_METRIC_PROCEDURES
  call debug_timer(stop_timer,"Reducing local and total volumes")
#endif
  !
  ! Close the files associated with output_metrics
  !
  if (output_metrics) then
#ifdef TIME_METRIC_PROCEDURES
    call debug_timer(start_accumulation,"output metrics",acc=9)
#endif
    do i = ibeg,iend
      close (iom(i),iostat=ierr,iomsg=error_message)
      call io_error(pname,fname(i),2,__LINE__,__FILE__,ierr,error_message)
    end do
#ifdef TIME_METRIC_PROCEDURES
    call debug_timer(stop_accumulation,"output metrics",acc=9)
#endif
  end if
  !
#ifdef TIME_METRIC_PROCEDURES
  call debug_timer(report_accumulated_time,"output metrics",acc=9)
#endif
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  21 format (/,60("-"),/," Cell ",i0,/,60("-"))
  22 format (5x,"Solution Point ",i4)
  23 format (3(2x,"d(",a,",",a,",",a,")/d",a,"  :  ",3es14.6,:,/))
  !
end subroutine old_get_metrics
!
!###############################################################################
!
subroutine old_compute_face_normals
  !
  !.. Use Statements ..
  use order_mod,         only : maxSP,maxFP
  use order_mod,         only : geom_solpts
  use geovar,            only : ncell,cell,fp
  use geovar,            only : nfbnd,nface,face,n_flxpts
  use geovar,            only : nr,xyz,bface
  use interpolation_mod, only : interp
  use derivatives_mod,   only : deriv
  !
  !.. Local Scalars ..
  integer :: i,j,m,n,nc,nf,ncf,f0,n1,n2,np
  integer :: kf,lf,l,k,nfp,ndim,lcfp,rcfp
  integer :: left_cell,left_geom,left_order
  integer :: right_cell,right_geom,right_order
  integer :: face_geom,face_order
  integer :: inrm,ierr,npcount,nfpts
  real(wp) :: a,val
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
  character(len=*), parameter :: pname = "old_compute_face_normals"
  character(len=*), parameter :: fname = "face_normals.tec.dat"
 !logical(lk), parameter :: output_face_normals = true
  logical(lk), parameter :: output_face_normals = fals
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
!#ifdef TIME_METRIC_PROCEDURES
!  call debug_timer(reset_accumulated_time,"output face normals",acc=1)
!  call debug_timer(reset_accumulated_time,"setup and transposing xyz",acc=2)
!  call debug_timer(reset_accumulated_time,"computing xyz_face",acc=3)
!  call debug_timer(reset_accumulated_time,"Filling individual geofa(nf)",acc=4)
!  call debug_timer(reset_accumulated_time,"Computing dvec from xyz_face",acc=5)
!  call debug_timer(reset_accumulated_time, &
!                   "Computing dvec cross products",acc=6)
!  call debug_timer(reset_accumulated_time,"Chopping dvec into geofa",acc=7)
!#endif
  !
!#ifdef TIME_METRIC_PROCEDURES
!  call debug_timer(start_timer,"Allocating geofa and its component arrays")
!#endif
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
!#ifdef TIME_METRIC_PROCEDURES
!  call debug_timer(stop_timer,"Allocating geofa and its component arrays")
!#endif
  !
  if (output_face_normals) then
!#ifdef TIME_METRIC_PROCEDURES
!    call debug_timer(start_accumulation,"output face normals",acc=1)
!#endif
    !
    npcount = 0
    nfpts = 0
    do nf = 1,nface
      nfpts = nfpts + size(geofa(nf)%v,dim=2)
    end do
    !
    if (nr > 1) then
      open (newunit=inrm,file=fname,status="replace",action="write", &
                         form="formatted",position="rewind",iostat=ierr, &
                         iomsg=error_message)
      call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
    end if
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
!#ifdef TIME_METRIC_PROCEDURES
!    call debug_timer(stop_accumulation,"output face normals",acc=1)
!#endif
  end if
  !
!#ifdef TIME_METRIC_PROCEDURES
!  call debug_timer(start_timer,"Computing geofa")
!#endif
  cell_loop: do nc = 1,ncell
    !
!#ifdef TIME_METRIC_PROCEDURES
!    call debug_timer(start_accumulation,"setup and transposing xyz",acc=2)
!#endif
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
!#ifdef TIME_METRIC_PROCEDURES
!    call debug_timer(stop_accumulation,"setup and transposing xyz",acc=2)
!#endif
    !
    ! Loop through the faces of the current cell
    ! to get the normals at each face point
    !
    cell_face_loop: do ncf = 1,size(cell(nc)%face)
      !
      if (cell(nc)%face(ncf)%side == 2) cycle cell_face_loop
!#ifdef TIME_METRIC_PROCEDURES
!      call debug_timer(start_accumulation,"computing xyz_face",acc=3)
!#endif
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
          lf = cell(nc)%face(ncf)%fp_offset + fp%geord(face_geom,face_order)% &
                                            rot(cell(nc)%face(ncf)%rotation)% &
                                              pts(kf)
          !
          val = zero
          do i = interp(left_geom,left_order)%toFace(face_order)%row_ptr(lf)+1,&
                 interp(left_geom,left_order)%toFace(face_order)%row_ptr(lf+1)
            !
            j = interp(left_geom,left_order)%toFace(face_order)%row_idx(i)
            a = interp(left_geom,left_order)%toFace(face_order)%row_val(i)
            !
            val = val + a * txyz(j,l)
            !
          end do
          !
          xyz_face(kf,l) = val
          !
        end do
        !
      end do
!#ifdef TIME_METRIC_PROCEDURES
!      call debug_timer(stop_accumulation,"computing xyz_face",acc=3)
!#endif
      !
      ! If the current face is an edge, initialize the k components of the
      ! two vectors the dvec array to one and everything else to zero
      !
!#ifdef TIME_METRIC_PROCEDURES
!      call debug_timer(start_accumulation,"Filling individual geofa(nf)",acc=4)
!#endif
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
!#ifdef TIME_METRIC_PROCEDURES
!        call debug_timer(start_accumulation, &
!                         "Computing dvec from xyz_face",acc=5)
!#endif
        do n = 1,ndim
          !
          ! If the current cell face is a x-min, y-max, or z-min face, we need
          ! to reverse the order of the derivatives so that the cross_product
          ! computed below results in a outward facing normal for the left cell
          ! on the face.
          !
          m = n
          if (face_geom /= Geom_Edge) then
            if (any(ncf == [Xmin3D_Face,Ymax3D_Face,Zmin3D_Face])) then
              m = ndim+1 - n
            end if
          end if
          !
          do l = 1,nr
            val = zero
            do i = deriv(face_geom,face_order)%d(n)%row_ptr(k)+1, &
                   deriv(face_geom,face_order)%d(n)%row_ptr(k+1)
              !
              j = deriv(face_geom,face_order)%d(n)%row_idx(i)
              a = deriv(face_geom,face_order)%d(n)%row_val(i)
              !
              val = val + a * xyz_face(j,l)
              !
            end do
            dvec(l,m) = val
          end do
          !
        end do
!#ifdef TIME_METRIC_PROCEDURES
!        call debug_timer(stop_accumulation, &
!                         "Computing dvec from xyz_face",acc=5)
!#endif
        !
        ! Finally, the length/area normal to the current face at the current
        ! flux point is just the cross product of the two vectors in dvec
        !
!#ifdef TIME_METRIC_PROCEDURES
!        call debug_timer(start_accumulation, &
!                         "Computing dvec cross products",acc=6)
!#endif
        dvec(:,3) = cross_product( dvec(:,1) , dvec(:,2) )
        !
        if (face_geom == Geom_Node) then
          dvec(1:nr,3) = one
         !dvec(1:nr,3) = merge(one,-one,ncf==1)
        end if
!#ifdef TIME_METRIC_PROCEDURES
!        call debug_timer(stop_accumulation, &
!                         "Computing dvec cross products",acc=6)
!#endif
        !
        ! Zero/chop all the normals less than 10*roundoff
        !
!#ifdef TIME_METRIC_PROCEDURES
!        call debug_timer(start_accumulation,"Chopping dvec into geofa",acc=7)
!#endif
        geofa(nf)%v(1:nr,k) = chop( real(dvec(1:nr,3),kind=qp)/qten ) * ten
!#ifdef TIME_METRIC_PROCEDURES
!        call debug_timer(stop_accumulation,"Chopping dvec into geofa",acc=7)
!#endif
        !
       !geofa(nf)%area = geofa(nf)%area + norm2(geofa(nf)%v(1:nr,k)) * &
       !                                  std_elem(face_geom,face_order)%wts(k)
        !
        if (output_face_normals) then
!#ifdef TIME_METRIC_PROCEDURES
!          call debug_timer(start_accumulation,"output face normals",acc=1)
!#endif
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
!#ifdef TIME_METRIC_PROCEDURES
!          call debug_timer(stop_accumulation,"output face normals",acc=1)
!#endif
        end if
        !
      end do
!#ifdef TIME_METRIC_PROCEDURES
!      call debug_timer(stop_accumulation,"Filling individual geofa(nf)",acc=4)
!#endif
      !
    end do cell_face_loop
    !
  end do cell_loop
!#ifdef TIME_METRIC_PROCEDURES
!  call debug_timer(stop_timer,"Computing geofa")
!#endif
  !
  ! Copy geofa into facenrm
  !
!#ifdef TIME_METRIC_PROCEDURES
!  call debug_timer(start_timer,"Migrating geofa(:)%v(:,:) to facenrm(:,:)")
!#endif
  allocate ( facenrm(1:nr,1:n_flxpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"facenrm",1,__LINE__,__FILE__,ierr,error_message)
  !$ACC UPDATE DEVICE(facenrm)
  !
  n = 0
  do nf = 1,size(geofa)
    do k = 1,size(geofa(nf)%v,dim=2)
      n = n + 1
      facenrm(1:nr,n) = geofa(nf)%v(1:nr,k)
    end do
  end do
  !$ACC UPDATE DEVICE(facenrm)
  !
  if (allocated(geofa)) then
    deallocate ( geofa , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"geofa",2,__LINE__,__FILE__,ierr,error_message)
  end if
!#ifdef TIME_METRIC_PROCEDURES
!  call debug_timer(stop_timer,"Migrating geofa(:)%v(:,:) to facenrm(:,:)")
!#endif
  !
  if (output_face_normals .and. nr == 1) then
!#ifdef TIME_METRIC_PROCEDURES
!    call debug_timer(start_accumulation,"output face normals",acc=1)
!#endif
    !
    write (array_name,'("face_normals.1D.CPU-",i2.2,".dat")') mypnum
    open (newunit=inrm,file=trim(adjustl(array_name)),status="replace", &
                       action="write",form="formatted",position="rewind", &
                       iostat=ierr,iomsg=error_message)
    call io_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    do nf = 1,nface
      !
      f0 = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      !
      write (inrm,30) nf
      do kf = 1,nfp
        write (inrm,31) f0+kf,kf,(facenrm(l,f0+kf), l=1,nr)
      end do
      !
    end do
    write (inrm,30)
    !
    close (inrm,iostat=ierr,iomsg=error_message)
    call io_error(pname,array_name,2,__LINE__,__FILE__,ierr,error_message)
    !
!#ifdef TIME_METRIC_PROCEDURES
!    call debug_timer(stop_accumulation,"output face normals",acc=1)
!#endif
  end if
  !
  30 format (80("="),:,/,2x,"Face ",i0,/,80("-"))
  31 format (4x," fp = ",i3,"; kf = ",i2,"; nrm = (",3(es14.6,")",:,tl1,","))
  !
  ! Create the interpolated jacobians at the face points
  !
  allocate ( facejac(1:n_flxpts,1:2) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"facejac",1,__LINE__,__FILE__,ierr,error_message)
  !$ACC UPDATE DEVICE(facejac)
  !
  ! Interpolate the jacobians from host cells
  ! to the face points of boundary faces
  !
!#ifdef TIME_METRIC_PROCEDURES
!  call debug_timer(start_timer,"Computing boundary face jacobians")
!#endif
  do nf = 1,nfbnd
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    left_cell  = face(nf)%left%cell
    left_geom  = cell(left_cell)%geom
    left_order = cell(left_cell)%order
    !
    do kf = 1,nfp
      !
      ! Interpolate the jacobian to the flux point for the left cell
      !
      lcfp = face(nf)%left%fp_offset + fp%geord(face_geom,face_order)% &
                                          rot(face(nf)%left%rotation)% &
                                       pts(kf)
      !
      val = zero
      do i = interp(left_geom,left_order)%toFace(face_order)%row_ptr(lcfp)+1, &
             interp(left_geom,left_order)%toFace(face_order)%row_ptr(lcfp+1)
        !
        j = interp(left_geom,left_order)%toFace(face_order)%row_idx(i)
        a = interp(left_geom,left_order)%toFace(face_order)%row_val(i)
        !
        val = val + a * metrics_dt(left_cell)%jac(j)
        !
      end do
      facejac(f0+kf,1) = val
      !
    end do
    !
  end do
!#ifdef TIME_METRIC_PROCEDURES
!  call debug_timer(stop_timer,"Computing boundary face jacobians")
!#endif
  !
  ! Now do the interior faces
  !
!#ifdef TIME_METRIC_PROCEDURES
!  call debug_timer(start_timer,"Computing interior face jacobians")
!#endif
  do nf = nfbnd+1,nface
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    left_cell  = face(nf)%left%cell
    right_cell = face(nf)%right%cell
    !
    left_geom  = cell(left_cell)%geom
    left_order = cell(left_cell)%order
    !
    right_geom  = cell(right_cell)%geom
    right_order = cell(right_cell)%order
    !
    do kf = 1,nfp
      !
      ! Interpolate the jacobian to the flux point for the left cell
      !
      lcfp = face(nf)%left%fp_offset + fp%geord(face_geom,face_order)% &
                                          rot(face(nf)%left%rotation)% &
                                       pts(kf)
      !
      val = zero
      do i = interp(left_geom,left_order)%toFace(face_order)%row_ptr(lcfp)+1, &
             interp(left_geom,left_order)%toFace(face_order)%row_ptr(lcfp+1)
        !
        j = interp(left_geom,left_order)%toFace(face_order)%row_idx(i)
        a = interp(left_geom,left_order)%toFace(face_order)%row_val(i)
        !
        val = val + a * metrics_dt(left_cell)%jac(j)
        !
      end do
      facejac(f0+kf,1) = val
      !
      ! Interpolate the jacobian to the flux point for the right cell
      !
      rcfp = face(nf)%right%fp_offset + fp%geord(face_geom,face_order)% &
                                          rot(face(nf)%right%rotation)% &
                                        pts(kf)
      !
      val = zero
      do i = interp(right_geom,right_order)%toFace(face_order)%row_ptr(rcfp)+1,&
             interp(right_geom,right_order)%toFace(face_order)%row_ptr(rcfp+1)
        !
        j = interp(right_geom,right_order)%toFace(face_order)%row_idx(i)
        a = interp(right_geom,right_order)%toFace(face_order)%row_val(i)
        !
        val = val + a * metrics_dt(right_cell)%jac(j)
        !
      end do
      facejac(f0+kf,2) = val
      !
    end do
    !
  end do
!#ifdef TIME_METRIC_PROCEDURES
!  call debug_timer(stop_timer,"Computing interior face jacobians")
!#endif
  !$ACC UPDATE DEVICE(facejac)
  !
  !
  !
  if (output_face_normals .and. nr > 1) then
!#ifdef TIME_METRIC_PROCEDURES
!    call debug_timer(start_accumulation,"output face normals",acc=1)
!#endif
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
!#ifdef TIME_METRIC_PROCEDURES
!    call debug_timer(stop_accumulation,"output face normals",acc=1)
!#endif
  end if
  !
!#ifdef TIME_METRIC_PROCEDURES
!  call debug_timer(report_accumulated_time,"output face normals",acc=1)
!  call debug_timer(report_accumulated_time,"setup and transposing xyz",acc=2)
!  call debug_timer(report_accumulated_time,"computing xyz_face",acc=3)
!  call debug_timer(report_accumulated_time, &
!                   "Filling individual geofa(nf)",acc=4)
!  call debug_timer(report_accumulated_time, &
!                   "Computing dvec from xyz_face",acc=5)
!  call debug_timer(report_accumulated_time, &
!                   "Computing dvec cross products",acc=6)
!  call debug_timer(report_accumulated_time,"Chopping dvec into geofa",acc=7)
!#endif
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
end subroutine old_compute_face_normals
!
!###############################################################################
!
subroutine compute_face_normals
  !
  !.. Use Statements ..
  use order_mod,         only : maxSP,maxFP
  use order_mod,         only : geom_solpts
  use geovar,            only : ncell,cell,fp
  use geovar,            only : nfbnd,nface,face,n_flxpts
  use geovar,            only : nr,xyz,bface
  use interpolation_mod, only : interp
  use derivatives_mod,   only : deriv
  !
  !.. Local Scalars ..
  integer :: i,j,m,n,nc,nf,ncf,f0,n1,n2,np
  integer :: kf,lf,l,k,nfp,ndim,lcfp,rcfp
  integer :: left_cell,left_geom,left_order
  integer :: right_cell,right_geom,right_order
  integer :: face_geom,face_order
  integer :: inrm,ierr,npcount,nfpts
  real(wp) :: a,val
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
#ifdef TIME_METRIC_PROCEDURES
  call debug_timer(reset_accumulated_time,"computing xyz_face",acc=1)
  call debug_timer(reset_accumulated_time, &
                   "Filling individual facenrm(:,k)",acc=2)
  call debug_timer(reset_accumulated_time,"Computing d(xyz)/di",acc=3)
  call debug_timer(reset_accumulated_time,"Computing d(xyz)/dj",acc=4)
  call debug_timer(reset_accumulated_time,"Computing cross product",acc=5)
  call debug_timer(reset_accumulated_time,"Chopping into facenrm",acc=6)
#endif
  !
  ! Allocate the facenrm
  !
  allocate ( facenrm(1:nr,1:n_flxpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"facenrm",1,__LINE__,__FILE__,ierr,error_message)
  !$ACC UPDATE DEVICE(facenrm)
  !
  ! Allocate face_rrf_vel
  !
  allocate ( face_rrf_vel(1:nr,1:n_flxpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"face_rrf_vel",1,__LINE__,__FILE__,ierr,error_message)
  !$ACC UPDATE DEVICE(face_rrf_vel)
  !
  ! Compute the face normals
  !
  if (nr == 3) then
    call compute_face_normals_3d(facenrm,face_rrf_vel)
  else if (nr == 2) then
    call compute_face_normals_2d(facenrm)
  else if (nr == 1) then
    call compute_face_normals_1d(facenrm)
  else
    write (error_message,1) nr
    call stop_gfr(stop_mpi,pname,__LINE__,__FILE__,error_message)
  end if
  !$ACC UPDATE DEVICE(facenrm,face_rrf_vel)
  !
  ! Create the interpolated jacobians at the face points
  !
  allocate ( facejac(1:n_flxpts,1:2) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"facejac",1,__LINE__,__FILE__,ierr,error_message)
  !$ACC UPDATE DEVICE(facejac)
  !
  ! Interpolate the cell metric jacobians to the face points
  !
  call compute_face_jacobians(facejac)
  !$ACC UPDATE DEVICE(facejac)
  !
#ifdef TIME_METRIC_PROCEDURES
  call debug_timer(report_accumulated_time,"computing xyz_face",acc=1)
  call debug_timer(report_accumulated_time, &
                   "Filling individual facenrm(:,k)",acc=2)
  call debug_timer(report_accumulated_time,"Computing d(xyz)/di",acc=3)
  call debug_timer(report_accumulated_time,"Computing d(xyz)/dj",acc=4)
  call debug_timer(report_accumulated_time,"Computing cross product",acc=5)
  call debug_timer(report_accumulated_time,"Chopping into facenrm",acc=6)
#endif
  !
  call debug_timer(leaving_procedure,pname)
  !
  ! Format Statements
  !
  1  format ("ERROR in subroutine compute_face_normals! The value of nr ", &
             "should be in the range of [1-3], but nr = ",i0,"!")
  !
end subroutine compute_face_normals
!
!###############################################################################
!
subroutine compute_face_normals_3d(fnrm,rrf_fvel)
  !
  !.. Use Statements ..
  use order_mod,         only : maxFP
  use order_mod,         only : geom_solpts
  use geovar,            only : nfbnd,ncell,cell,fp,face
  use geovar,            only : xyz,xyz_nodes
  use geovar,            only : nodes_of_cell,nodes_of_cell_ptr
  use interpolation_mod, only : interp
  use derivatives_mod,   only : deriv
  !
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  use omp_lib, only : omp_get_thread_num
  use omp_lib, only : omp_get_num_threads
#endif
#endif
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: fnrm(:,:)
  real(wp), intent(inout) :: rrf_fvel(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,n,nc,nf,ncf,n0
  integer  :: kf,lf,k,f0,nfp
  integer  :: this_geom,this_order
  integer  :: face_geom,face_order
  integer  :: tc_offset,tc_rot
  real(wp) :: a,facenrm_coeff
  real(wp) :: dxdi,dydi,dzdi
  real(wp) :: dxdj,dydj,dzdj
  real(wp) :: xval,yval,zval
  real(wp) :: xnrm,ynrm,znrm
  real(wp) :: x,y,z,r,theta,angle
  real(qp) :: xqp,yqp,zqp
  !
  !.. Local Arrays ..
  real(wp) :: xyz_face(1:maxFP,1:3)
  !
  integer :: i1,i2
  integer, allocatable :: cell_nodes(:)
  integer, allocatable :: face_nodes(:)
  real(wp), allocatable :: cell_xyz(:,:)
  real(wp), allocatable :: face_xyz(:,:)
  real(wp) :: cell_vec(1:3)
  real(wp) :: cell_pt(1:3)
  real(wp) :: face_pt(1:3)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_face_normals_3d"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
!#undef TIME_METRIC_PROCEDURES
  !
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  !$OMP PARALLEL &
  !$OMP&         DEFAULT(SHARED) &
  !$OMP&         PRIVATE(nc,this_geom,this_order,n0,ncf,nf, &
  !$OMP&                 face_geom,face_order,f0,nfp,tc_offset,tc_rot, &
  !$OMP&                 kf,lf,xval,yval,zval,i,k,n,a,xyz_face,tid,  &
  !$OMP&                 dxdi,dxdj,dydi,dydj,dzdi,dzdj, &
  !$OMP&                 xnrm,ynrm,znrm,xqp,yqp,zqp)
  !
  tid = omp_get_thread_num()
  nthreads = omp_get_num_threads()
  !
  if (mypnum == glb_root) then
    if (tid == glb_root) then
      write (iout,thread_report) pname,nthreads
    end if
  end if
  !
  chunk_size = ceiling( real(ncell,kind=wp) / real(nthreads,kind=wp) )
  chunk_size = chunk_size / 10
  !
  !$OMP DO SCHEDULE(RUNTIME)
 !!$OMP DO SCHEDULE(GUIDED,chunk_size)
 !!$OMP DO SCHEDULE(STATIC,chunk_size)
#endif
#endif
  cell_loop: do nc = 1,ncell
    !
    this_geom  = cell(nc)%geom
    this_order = cell(nc)%order
    !
    n0 = cell(nc)%beg_sp - 1
    !
    i1 = nodes_of_cell_ptr(nc)+1
    i2 = nodes_of_cell_ptr(nc+1)
    !
    cell_nodes = nodes_of_cell(i1:i2) ! F2003 AUTO-REALLOC
    cell_xyz = xyz_nodes(1:3,cell_nodes(:)) ! F2003 AUTO-REALLOC
    !
    ! Loop through the faces of the current cell
    ! to get the normals at each face point
    !
    cell_face_loop: do ncf = 1,size(cell(nc)%face)
      !
      if (cell(nc)%face(ncf)%side == 2) cycle cell_face_loop
#ifdef TIME_METRIC_PROCEDURES
      call debug_timer(start_accumulation,"computing xyz_face",acc=1)
#endif
      !
      nf = cell(nc)%face(ncf)%idx
      !
      face_geom  = cell(nc)%face(ncf)%geom
      face_order = cell(nc)%face(ncf)%order
      !
     !nfp = geom_solpts(face_geom,face_order)
      f0 = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      !
      tc_offset = cell(nc)%face(ncf)%fp_offset
      tc_rot = cell(nc)%face(ncf)%rotation
      !
      if (face_geom == Geom_Edge) then
        face_nodes = face(nf)%nodes(1:2) ! F2003 AUTO-REALLOC
      else if (face_geom == Geom_Tria) then
        face_nodes = face(nf)%nodes(1:3) ! F2003 AUTO-REALLOC
      else
        face_nodes = face(nf)%nodes(1:4) ! F2003 AUTO-REALLOC
      end if
      face_xyz = xyz_nodes(1:3,face_nodes(:)) ! F2003 AUTO-REALLOC
      !
      call get_cell_vec_sub(cell_nodes,face_nodes,cell_xyz,face_xyz, &
                            cell_pt,face_pt,cell_vec)
     !cell_vec(1:3) = get_cell_vec(cell_nodes,face_nodes,cell_xyz,face_xyz)
      !
      ! Interpolate the coordinates of the solution points to get the
      ! coordinates of the flux points on the current face.
      !
      do kf = 1,nfp
        !
        lf = tc_offset + fp%geord(face_geom,face_order)%rot(tc_rot)%pts(kf)
        !
        xval = zero
        yval = zero
        zval = zero
        !
        do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1,&
               interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
          !
          n = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
          a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
          !
          xval = xval + a * xyz(1,n0+n)
          yval = yval + a * xyz(2,n0+n)
          zval = zval + a * xyz(3,n0+n)
          !
        end do
        !
        xyz_face(kf,1) = xval
        xyz_face(kf,2) = yval
        xyz_face(kf,3) = zval
        !
      end do
#ifdef TIME_METRIC_PROCEDURES
      call debug_timer(stop_accumulation,"computing xyz_face",acc=1)
#endif
#ifdef TIME_METRIC_PROCEDURES
      call debug_timer(start_accumulation, &
                       "Filling individual facenrm(:,k)",acc=2)
#endif
      !
      ! We need to create a coefficient to reverse the direction of the
      ! face normal depending on which face of the current cell. This
      ! coefficient will have the default value of +1, but will be changed
      ! to -1 for the faces that need the face normal direction reversed.
      !
      facenrm_coeff = one
      !
      if (this_geom == Geom_Hexa) then
        !
        ! If the current cell is a hexahedra, and the face is a x-min, y-max,
        ! or z-min face, we need to reverse the order of the derivatives so
        ! that the cross_product computed below results in a outward facing
        ! normal for the left cell on the face. This means we compute
        ! [d(xyz)dj]x[d(xyz)di].
        !
        if (any(ncf == [Xmin3D_Face,Ymax3D_Face,Zmin3D_Face])) then
          facenrm_coeff = -one
        end if
        !
      else if (this_geom == Geom_Tetr) then
        !
        ! Faces 3 and 4 of tetrahedras need the face normal direction reversed
        !
        if (any(ncf == [3,4])) then
          facenrm_coeff = -one
        end if
        !
      else if (this_geom == Geom_Pris) then
        !
        ! Faces 1 and 5 of prisms need the face normal direction reversed
        !
        if (any(ncf == [1,5])) then
          facenrm_coeff = -one
        end if
        !
      else if (this_geom == Geom_Pyra) then
        !
        ! Faces 3, 4, and 5 of pyramids need the face normal direction reversed
        !
        if (any(ncf == [3,4,5])) then
          facenrm_coeff = -one
        end if
        !
      end if
      !
      ! Now compute the derivative of the flux point coordinates
      !
      do kf = 1,nfp
        !
        dxdi = zero ; dxdj = zero
        dydi = zero ; dydj = zero
        dzdi = zero ; dzdj = zero
        !
        ! Xi/I derivative of the face
        !
#ifdef TIME_METRIC_PROCEDURES
        call debug_timer(start_accumulation,"Computing d(xyz)/di",acc=3)
#endif
        do i = deriv(face_geom,face_order)%d(1)%row_ptr(kf)+1, &
               deriv(face_geom,face_order)%d(1)%row_ptr(kf+1)
          !
          k = deriv(face_geom,face_order)%d(1)%row_idx(i)
          a = deriv(face_geom,face_order)%d(1)%row_val(i)
          !
          dxdi = dxdi + a * xyz_face(k,1)
          dydi = dydi + a * xyz_face(k,2)
          dzdi = dzdi + a * xyz_face(k,3)
          !
        end do
#ifdef TIME_METRIC_PROCEDURES
        call debug_timer(stop_accumulation,"Computing d(xyz)/di",acc=3)
#endif
        !
        ! Eta/J derivative of the face
        !
#ifdef TIME_METRIC_PROCEDURES
        call debug_timer(start_accumulation,"Computing d(xyz)/dj",acc=4)
#endif
        do i = deriv(face_geom,face_order)%d(2)%row_ptr(kf)+1, &
               deriv(face_geom,face_order)%d(2)%row_ptr(kf+1)
          !
          k = deriv(face_geom,face_order)%d(2)%row_idx(i)
          a = deriv(face_geom,face_order)%d(2)%row_val(i)
          !
          dxdj = dxdj + a * xyz_face(k,1)
          dydj = dydj + a * xyz_face(k,2)
          dzdj = dzdj + a * xyz_face(k,3)
          !
        end do
#ifdef TIME_METRIC_PROCEDURES
        call debug_timer(stop_accumulation,"Computing d(xyz)/dj",acc=4)
#endif
        !
        ! Compute [d(xyz)di] x [d(xyz)dj]
        ! NOTE: If facenrm_coeff = -1, this will compute [d(xyz)dj] x [d(xyz)di]
        !
#ifdef TIME_METRIC_PROCEDURES
        call debug_timer(start_accumulation,"Computing cross product",acc=5)
#endif
        xnrm = facenrm_coeff * (dydi*dzdj - dzdi*dydj)
        ynrm = facenrm_coeff * (dzdi*dxdj - dxdi*dzdj)
        znrm = facenrm_coeff * (dxdi*dydj - dydi*dxdj)
#ifdef TIME_METRIC_PROCEDURES
        call debug_timer(stop_accumulation,"Computing cross product",acc=5)
#endif
        !
        ! Zero/chop all the normals less than 10*roundoff
        ! NOTE: Using the chop() function here was really slow.
        !       Manually inlining the function here resulted in
        !       over two orders of magnitude speedup for this
        !       subroutine.
        !
#ifdef TIME_METRIC_PROCEDURES
        call debug_timer(start_accumulation,"Chopping into facenrm",acc=6)
#endif
        xqp = real(xnrm,kind=qp)
        yqp = real(ynrm,kind=qp)
        zqp = real(znrm,kind=qp)
        !
        xqp = merge( xqp , qzero , abs(xqp/qten)>real(rndoff,kind=qp) )
        yqp = merge( yqp , qzero , abs(yqp/qten)>real(rndoff,kind=qp) )
        zqp = merge( zqp , qzero , abs(zqp/qten)>real(rndoff,kind=qp) )
        !
        fnrm(1,f0+kf) = real(xqp,kind=wp)
        fnrm(2,f0+kf) = real(yqp,kind=wp)
        fnrm(3,f0+kf) = real(zqp,kind=wp)
#ifdef TIME_METRIC_PROCEDURES
        call debug_timer(stop_accumulation,"Chopping into facenrm",acc=6)
#endif
        !
        if (dot(fnrm(1:3,f0+kf),cell_vec(1:3)) > zero) then
       !if (this_geom == Geom_Pyra) then
          ! If any of the output dot products are positive, that face
          ! of the cell (ncf) needs to have logic to set facenrm_coeff
          ! to -one. The logic is correctly in place above for faces
          ! of hexes, tets, pyramids, and prisms.
          x = xyz_face(kf,1)
          y = xyz_face(kf,2)
          z = xyz_face(kf,3)
          r = sqrt(y*y + z*z)
          theta = atan2(z,y)
          angle = (theta + merge(two*pi,zero,theta<zero))*radians2degrees
          write (9101,95) nc,this_geom,f0+kf,kf,ncf, &
                          dot(fnrm(1:3,f0+kf),cell_vec(1:3)), &
                          x,y,z,r,theta,angle, &
                          (fnrm(i,f0+kf),i=1,3),(cell_vec(i),i=1,3), &
                          (cell_pt(i),i=1,3),(face_pt(i),i=1,3)
        end if
        95 format (i6,i3,i10,i4,i8,f15.10,6(4x,3f11.6))
        !
        ! Compute the angular velocity vector of the
        ! rotating reference frame at this flux point
        !
        rrf_fvel(1:3,f0+kf) = rrf_face_velocity(xyz_face(kf,1:3))
        !
      end do
#ifdef TIME_METRIC_PROCEDURES
      call debug_timer(stop_accumulation, &
                       "Filling individual facenrm(:,k)",acc=2)
#endif
!#define TIME_METRIC_PROCEDURES
      !
    end do cell_face_loop
    !
  end do cell_loop
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  !$OMP END DO
  !$OMP END PARALLEL
#endif
#endif
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_face_normals_3d
!
!###############################################################################
!
subroutine get_cell_vec_sub(cell_nodes,face_nodes,cell_xyz,face_xyz, &
                            cell_pt,face_pt,cell_vec)
  !
  integer,     intent(in) :: cell_nodes(:)
  integer,     intent(in) :: face_nodes(:)
  real(wp),    intent(in) :: cell_xyz(:,:)
  real(wp),    intent(in) :: face_xyz(:,:)
  real(wp), intent(inout) :: cell_pt(1:size(cell_xyz,dim=1))
  real(wp), intent(inout) :: face_pt(1:size(cell_xyz,dim=1))
  real(wp), intent(inout) :: cell_vec(1:size(cell_xyz,dim=1))
  !
  integer :: n,l,nfp,ncp,nl,m
  !
  real(wp) :: face_center(1:size(cell_xyz,dim=1))
  !
continue
  !
  cell_vec(:) = -huge(zero)
  !
  nl = size(cell_xyz,dim=1)
  nfp = size(face_nodes)
  ncp = size(cell_nodes)
  !
  ! Compute the coordinates for the center of the current face by averaging
  ! the coordinates of all the face nodes.
  !
  do l = 1,nl
    face_center(l) = sum(face_xyz(l,:)) / real(nfp,kind=wp)
  end do
  !
  ! Loop through the cell nodes to find a node that is not on the current face
  !
  m = 0
  cell_pt(:) = zero
  do n = 1,ncp
    !
    ! If this cell node matches any of the face nodes, we need to skip to
    ! the next cell node because we need to find all the cell nodes that are
    ! not on the current face to get a vector from the cell interior (direction)
    ! to the center of the current face.
    !
    if (any(cell_nodes(n) == face_nodes(:))) cycle
    !
    ! This cell node is not part of the face, so add it to the average of
    ! the nodes not on the face
    !
    m = m + 1
    !
    do l = 1,nl
      cell_pt(l) = cell_pt(l) + cell_xyz(l,n)
    end do
    !
  end do
  !
  ! Now compute the average location of all nodes not on the face using
  ! the counter m
  !
  do l = 1,nl
    cell_pt(l) = cell_pt(l)/real(m,kind=wp)
  end do
  !
  ! Finally, compute the vector from the averaged coordinates
  ! of the cell nodes not on the face, to the center of the face
  !
  do l = 1,nl
    cell_vec(l) = cell_pt(l) - face_center(l)
    face_pt(l) = face_center(l)
  end do
 !!
 !! Loop through the cell nodes to find a node that is not on the current face
 !!
 !do n = 1,ncp
 !  !
 !  ! If this cell node matches any of the face nodes, we need to skip to
 !  ! the next cell node because we need to a cell node that is not on the
 !  ! current face to get a vector from the cell interior (direction) to
 !  ! the center of the current face.
 !  !
 !  if (any(cell_nodes(n) == face_nodes(:))) cycle
 !  !
 !  ! This cell node is not on the current face, so compute the vector from
 !  ! this cell node to the center of the face, which is a vector from the
 !  ! interior (direction) of the cell to the face.
 !  !
 !  do l = 1,nl
 !    cell_vec(l) = cell_xyz(l,n) - face_center(l)
 !    cell_pt(l) = cell_xyz(l,n)
 !    face_pt(l) = face_center(l)
 !  end do
 !  !
 !  ! We only need to compute this vector from one of the cell nodes that is
 !  ! not on the face, so if we are here, the job is done and we can exit.
 !  !
 !  exit
 !  !
 !end do
  !
end subroutine get_cell_vec_sub
!
!###############################################################################
!
function get_cell_vec(cell_nodes,face_nodes,cell_xyz,face_xyz) &
               result(return_value)
  !
  integer,  intent(in) :: cell_nodes(:)
  integer,  intent(in) :: face_nodes(:)
  real(wp), intent(in) :: cell_xyz(:,:)
  real(wp), intent(in) :: face_xyz(:,:)
  !
  real(wp) :: return_value(1:size(cell_xyz,dim=1))
  !
  integer :: n,l,nfp,ncp,nl
  !
  real(wp) :: face_center(1:size(cell_xyz,dim=1))
  !
continue
  !
  return_value(:) = -huge(zero)
  !
  nl = size(cell_xyz,dim=1)
  nfp = size(face_nodes)
  ncp = size(cell_nodes)
  !
  ! Compute the coordinates for the center of the current face by averaging
  ! the coordinates of all the face nodes.
  !
  do l = 1,nl
    face_center(l) = sum(face_xyz(l,:)) / real(nfp,kind=wp)
  end do
  !
  ! Loop through the cell nodes to find a node that is not on the current face
  !
  do n = 1,ncp
    !
    ! If this cell node matches any of the face nodes, we need to skip to
    ! the next cell node because we need to a cell node that is not on the
    ! current face to get a vector from the cell interior (direction) to
    ! the center of the current face.
    !
    if (any(cell_nodes(n) == face_nodes(:))) cycle
    !
    ! This cell node is not on the current face, so compute the vector from
    ! this cell node to the center of the face, which is a vector from the
    ! interior (direction) of the cell to the face.
    !
    do l = 1,nl
      return_value(l) = cell_xyz(l,n) - face_center(l)
    end do
    !
    ! We only need to compute this vector from one of the cell nodes that is
    ! not on the face, so if we are here, the job is done and we can exit.
    !
    exit
    !
  end do
  !
end function get_cell_vec
!
!###############################################################################
!
subroutine compute_face_normals_2d(fnrm)
  !
  !.. Use Statements ..
  use order_mod,         only : maxFP
  use order_mod,         only : geom_solpts
  use geovar,            only : ncell,cell,fp,face
  use geovar,            only : xyz
  use interpolation_mod, only : interp
  use derivatives_mod,   only : deriv
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: fnrm(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,n,nc,nf,ncf,n0
  integer  :: kf,lf,k,f0,nfp
  integer  :: this_geom,this_order
  integer  :: face_geom,face_order
  integer  :: tc_offset,tc_rot
  real(wp) :: a
  real(wp) :: dxdi,dydi
  real(wp) :: xval,yval
  !
  !.. Local Arrays ..
  real(wp) :: xyz_face(1:maxFP,1:2)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_face_normals_2d"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  cell_loop: do nc = 1,ncell
    !
    this_geom  = cell(nc)%geom
    this_order = cell(nc)%order
    !
    n0 = cell(nc)%beg_sp - 1
    !
    ! Loop through the faces of the current cell
    ! to get the normals at each face point
    !
    cell_face_loop: do ncf = 1,size(cell(nc)%face)
      !
      if (cell(nc)%face(ncf)%side == 2) cycle cell_face_loop
#ifdef TIME_METRIC_PROCEDURES
      call debug_timer(start_accumulation,"computing xyz_face",acc=1)
#endif
      !
      nf = cell(nc)%face(ncf)%idx
      !
      face_geom  = cell(nc)%face(ncf)%geom
      face_order = cell(nc)%face(ncf)%order
      !
     !nfp = geom_solpts(face_geom,face_order)
      f0 = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      !
      tc_offset = cell(nc)%face(ncf)%fp_offset
      tc_rot = cell(nc)%face(ncf)%rotation
      !
      ! Interpolate the coordinates of the solution points to get the
      ! coordinates of the flux points on the current face.
      !
      do kf = 1,nfp
        !
        lf = tc_offset + fp%geord(face_geom,face_order)%rot(tc_rot)%pts(kf)
        !
        xval = zero
        yval = zero
        !
        do i = interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf)+1,&
               interp(this_geom,this_order)%toFace(face_order)%row_ptr(lf+1)
          !
          n = interp(this_geom,this_order)%toFace(face_order)%row_idx(i)
          a = interp(this_geom,this_order)%toFace(face_order)%row_val(i)
          !
          xval = xval + a * xyz(1,n0+n)
          yval = yval + a * xyz(2,n0+n)
          !
        end do
        !
        xyz_face(kf,1) = xval
        xyz_face(kf,2) = yval
        !
      end do
#ifdef TIME_METRIC_PROCEDURES
      call debug_timer(stop_accumulation,"computing xyz_face",acc=1)
#endif
#ifdef TIME_METRIC_PROCEDURES
      call debug_timer(start_accumulation, &
                       "Filling individual facenrm(:,k)",acc=2)
#endif
      !
      ! Now compute the derivative of the flux point coordinates
      !
      do kf = 1,nfp
        !
        ! Xi/I derivative of the face
        !
#ifdef TIME_METRIC_PROCEDURES
        call debug_timer(start_accumulation,"Computing d(xyz)/di",acc=3)
#endif
        !
        dxdi = zero
        dydi = zero
        !
        do i = deriv(face_geom,face_order)%d(1)%row_ptr(kf)+1, &
               deriv(face_geom,face_order)%d(1)%row_ptr(kf+1)
          !
          k = deriv(face_geom,face_order)%d(1)%row_idx(i)
          a = deriv(face_geom,face_order)%d(1)%row_val(i)
          !
          dxdi = dxdi + a * xyz_face(k,1)
          dydi = dydi + a * xyz_face(k,2)
          !
        end do
#ifdef TIME_METRIC_PROCEDURES
        call debug_timer(stop_accumulation,"Computing d(xyz)/di",acc=3)
#endif
        !
        ! NOTE: We dont need to do a cross product for 2D since the
        !       x-normal is +dy/di and the y-normal is -dx/di.
        !
        ! Zero/chop all the normals less than 10*roundoff
        !
#ifdef TIME_METRIC_PROCEDURES
        call debug_timer(start_accumulation,"Chopping into facenrm",acc=6)
#endif
        fnrm(1,f0+kf) = chop( real( dydi,kind=qp)/qten ) * ten
        fnrm(2,f0+kf) = chop( real(-dxdi,kind=qp)/qten ) * ten
#ifdef TIME_METRIC_PROCEDURES
        call debug_timer(stop_accumulation,"Chopping into facenrm",acc=6)
#endif
        !
      end do
#ifdef TIME_METRIC_PROCEDURES
      call debug_timer(stop_accumulation, &
                       "Filling individual facenrm(:,k)",acc=2)
#endif
      !
    end do cell_face_loop
    !
  end do cell_loop
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_face_normals_2d
!
!###############################################################################
!
subroutine compute_face_normals_1d(fnrm)
  !
  !.. Use Statements ..
  use geovar, only : ncell,cell,face
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: fnrm(:,:)
  !
  !.. Local Scalars ..
  integer :: nc,nf,ncf,f0
  integer :: kf,nfp
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_face_normals_1d"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! NOTE: This whole subroutine can probably be simplified to:
  !       fnrm(:,:) = one
  !
  cell_loop: do nc = 1,ncell
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
      f0 = face(nf)%beg_fp - 1
      nfp = face(nf)%end_fp - f0
      !
      ! Now compute the derivative of the flux point coordinates
      !
      do kf = 1,nfp
        !
        ! The face normal for 1D cells is always +1
        !
        fnrm(1,f0+kf) = one
        !
      end do
      !
    end do cell_face_loop
    !
  end do cell_loop
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_face_normals_1d
!
!###############################################################################
!
pure function rrf_face_velocity(xyz) result(return_value)
  !
  !.. Use Statements ..
  use ovar, only : aref
  use ovar, only : ref_frame_rotation_vector
  !
  !.. Formal Arguments ..
  real(wp), intent(in) :: xyz(:)
  !
  !.. Function Result ..
  real(wp) :: return_value(1:size(xyz))
  !
  !.. Local Scalars ..
  integer  :: n
  real(wp) :: x,x_r
  real(wp) :: y,y_r
  real(wp) :: z,z_r
  !
  !.. Local Arrays ..
  real(wp) :: av(1:3)
  !
continue
  !
  return_value(:) = zero
  !
  n = size(xyz)
  !
  y = zero
  z = zero
  !
  ! NOTE: ref_frame_rotation_vector is the rotation rate vector given in
  !       RPM. We are pre-multiplying this through by rpm_to_mps/aref to
  !       convert RPM to m/s and then nondimensionalize by the reference
  !       speed of sound. Technically, this is converting RPM to rad/s, but
  !       we will be taking the cross-product of this rotation rate vector
  !       with the position vector which will convert from rad/s to m/s.
  !
  x_r = ref_frame_rotation_vector(1)*rpm_to_mps/aref ! x rotation
  y_r = ref_frame_rotation_vector(2)*rpm_to_mps/aref ! y rotation
  z_r = ref_frame_rotation_vector(3)*rpm_to_mps/aref ! z rotation
  !
  x = xyz(1)
  if (n > 1) y = xyz(2)
  if (n > 2) z = xyz(3)
  !
  av(1) = y_r*z - z_r*y ! angular velocity in x direction
  av(2) = z_r*x - x_r*z ! angular velocity in y direction
  av(3) = x_r*y - y_r*x ! angular velocity in z direction
  !
  return_value(1:n) = av(1:n)
  !
end function rrf_face_velocity
!
!###############################################################################
!
subroutine compute_face_jacobians(fjac)
  !
  !.. Use Statements ..
  use geovar,            only : nfbnd,nface,face,cell,fp
  use interpolation_mod, only : interp
  !
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  use omp_lib, only : omp_get_thread_num
  use omp_lib, only : omp_get_num_threads
#endif
#endif
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: fjac(:,:)
  !
  !.. Local Scalars ..
  integer  :: i,n,nf,kf,f0,nfp
  integer  :: l0,lc_offset,lc_rot,kl
  integer  :: r0,rc_offset,rc_rot,kr
  integer  :: host_cell,host_geom,host_order
  integer  :: left_cell,left_geom,left_order
  integer  :: right_cell,right_geom,right_order
  integer  :: face_geom,face_order
  real(wp) :: a,val
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_face_jacobians"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
#ifdef TIME_METRIC_PROCEDURES
  call debug_timer(start_timer,"Computing boundary face jacobians")
#endif
  !
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  !$OMP PARALLEL &
  !$OMP&         DEFAULT(SHARED) &
  !$OMP&         PRIVATE(tid,nf,face_geom,face_order,f0,nfp, &
  !$OMP&                 host_cell,host_geom,host_order, &
  !$OMP&                 l0,lc_offset,lc_rot,kl, &
  !$OMP&                 kf,val,i,n,a)
  !
  tid = omp_get_thread_num()
  nthreads = omp_get_num_threads()
  !
  if (mypnum == glb_root) then
    if (tid == glb_root) then
      write (iout,thread_report) pname,nthreads
    end if
  end if
  !
  chunk_size = ceiling( real(nfbnd,kind=wp) / real(nthreads,kind=wp) )
  chunk_size = chunk_size / 10
  !
  !$OMP DO SCHEDULE(RUNTIME)
 !!$OMP DO SCHEDULE(GUIDED,chunk_size)
 !!$OMP DO SCHEDULE(STATIC,chunk_size)
#endif
#endif
  do nf = 1,nfbnd
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    host_cell  = face(nf)%left%cell
    !
    host_geom  = cell(host_cell)%geom
    host_order = cell(host_cell)%order
    !
    l0 = cell(host_cell)%beg_sp-1
    !
    lc_offset = face(nf)%left%fp_offset
    lc_rot = face(nf)%left%rotation
    !
    do kf = 1,nfp
      !
      ! Interpolate the jacobian to the flux point for the left cell
      !
      kl = lc_offset + fp%geord(face_geom,face_order)%rot(lc_rot)%pts(kf)
      !
      val = zero
      do i = interp(host_geom,host_order)%toFace(face_order)%row_ptr(kl)+1, &
             interp(host_geom,host_order)%toFace(face_order)%row_ptr(kl+1)
        !
        n = interp(host_geom,host_order)%toFace(face_order)%row_idx(i)
        a = interp(host_geom,host_order)%toFace(face_order)%row_val(i)
        !
        val = val + a * jac(l0+n)
        !
      end do
      fjac(f0+kf,1) = val
      !
    end do
    !
  end do
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  !$OMP END DO
  !$OMP END PARALLEL
#endif
#endif
#ifdef TIME_METRIC_PROCEDURES
  call debug_timer(stop_timer,"Computing boundary face jacobians")
#endif
  !
  ! Now do the interior faces
  !
#ifdef TIME_METRIC_PROCEDURES
  call debug_timer(start_timer,"Computing interior face jacobians")
#endif
  !
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  !$OMP PARALLEL &
  !$OMP&         DEFAULT(SHARED) &
  !$OMP&         PRIVATE(tid,nf,face_geom,face_order,f0,nfp, &
  !$OMP&                 left_cell,left_geom,left_order, &
  !$OMP&                 l0,lc_offset,lc_rot,kl, &
  !$OMP&                 right_cell,right_geom,right_order, &
  !$OMP&                 r0,rc_offset,rc_rot,kr, &
  !$OMP&                 kf,val,i,n,a)
  !
  tid = omp_get_thread_num()
  nthreads = omp_get_num_threads()
  !
  if (mypnum == glb_root) then
    if (tid == glb_root) then
      write (iout,thread_report) pname,nthreads
    end if
  end if
  !
  chunk_size = ceiling( real(nface-nfbnd,kind=wp) / real(nthreads,kind=wp) )
  chunk_size = chunk_size / 10
  !
  !$OMP DO SCHEDULE(RUNTIME)
 !!$OMP DO SCHEDULE(GUIDED,chunk_size)
 !!$OMP DO SCHEDULE(STATIC,chunk_size)
#endif
#endif
  do nf = nfbnd+1,nface
    !
    face_geom = face(nf)%geom
    face_order = face(nf)%order
    !
    f0 = face(nf)%beg_fp - 1
    nfp = face(nf)%end_fp - f0
    !
    left_cell  = face(nf)%left%cell
    right_cell = face(nf)%right%cell
    !
    left_geom  = cell(left_cell)%geom
    left_order = cell(left_cell)%order
    !
    right_geom  = cell(right_cell)%geom
    right_order = cell(right_cell)%order
    !
    l0 = cell(left_cell)%beg_sp-1
    r0 = cell(right_cell)%beg_sp-1
    !
    lc_offset = face(nf)%left%fp_offset
    lc_rot = face(nf)%left%rotation
    !
    rc_offset = face(nf)%right%fp_offset
    rc_rot = face(nf)%right%rotation
    !
    do kf = 1,nfp
      !
      ! Interpolate the jacobian to the flux point for the left cell
      !
      kl = lc_offset + fp%geord(face_geom,face_order)%rot(lc_rot)%pts(kf)
      !
      val = zero
      do i = interp(left_geom,left_order)%toFace(face_order)%row_ptr(kl)+1, &
             interp(left_geom,left_order)%toFace(face_order)%row_ptr(kl+1)
        !
        n = interp(left_geom,left_order)%toFace(face_order)%row_idx(i)
        a = interp(left_geom,left_order)%toFace(face_order)%row_val(i)
        !
        val = val + a * jac(l0+n)
        !
      end do
      fjac(f0+kf,1) = val
      !
      ! Interpolate the jacobian to the flux point for the right cell
      !
      kr = rc_offset + fp%geord(face_geom,face_order)%rot(rc_rot)%pts(kf)
      !
      val = zero
      do i = interp(right_geom,right_order)%toFace(face_order)%row_ptr(kr)+1,&
             interp(right_geom,right_order)%toFace(face_order)%row_ptr(kr+1)
        !
        n = interp(right_geom,right_order)%toFace(face_order)%row_idx(i)
        a = interp(right_geom,right_order)%toFace(face_order)%row_val(i)
        !
        val = val + a * jac(r0+n)
        !
      end do
      fjac(f0+kf,2) = val
      !
    end do
    !
  end do
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  !$OMP END DO
  !$OMP END PARALLEL
#endif
#endif
#ifdef TIME_METRIC_PROCEDURES
  call debug_timer(stop_timer,"Computing interior face jacobians")
#endif
  !
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_face_jacobians
!
!###############################################################################
!
subroutine compute_quadrature_face_normals
  !
  !.. Use Statements ..
  use order_mod,         only : maxSP,maxGP
  use order_mod,         only : geom_solpts
  use geovar,            only : ncell,cell,fp
  use geovar,            only : nface,face,n_flx_qdtpts
  use geovar,            only : nr,xyz
  use interpolation_mod, only : interp
  use derivatives_mod,   only : deriv
  !
  !.. Local Scalars ..
  integer :: i,j,m,n,nc,nf,ncf,n0,np
  integer :: kf,lf,l,k,nqp,ndim
  integer :: left_geom,left_order
  integer :: face_geom,qdtr_order
  integer :: inrm,ierr,npcount,nfqpts
  real(wp) :: a,val
  character(len=200) :: array_name
  !
  !.. Local Arrays ..
  real(wp) :: dvec(1:3,1:3)
  real(wp) :: txyz(1:maxSP,1:nr)
  real(wp) :: xyz_qdtr(1:maxGP,1:nr)
  !
  !.. Local Allocatable Arrays ..
  integer,          allocatable :: ipface(:,:)
  type(face_nrm_t), allocatable :: geofq(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_quadrature_face_normals"
  character(len=*), parameter :: fname = "face_qdtr_normals.tec.dat"
 !logical(lk), parameter :: output_face_normals = true
  logical(lk), parameter :: output_face_normals = fals
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  ! Allocate the outer part of the geofq array
  !
  allocate ( geofq(1:nface) , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"geofq",1,__LINE__,__FILE__,ierr,error_message)
  !
  ! Allocate the data components of the geofq array
  !
  do nf = 1,nface
    !
    face_geom = face(nf)%geom
    qdtr_order = face(nf)%qdtr_order
    !
    nfqpts = geom_solpts(face_geom,qdtr_order)
    !
    allocate ( geofq(nf)%v(1:nr,1:nfqpts) , source=zero , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) "geofq",nf
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr, &
                     error_message,skip_alloc_pause)
    !
  end do
  !
  if (output_face_normals) then
    !
    npcount = 0
    nfqpts = 0
    do nf = 1,nface
      nfqpts = nfqpts + size(geofq(nf)%v,dim=2)
    end do
    !
    open (newunit=inrm,file=fname,status="replace",action="write", &
                       form="formatted",position="rewind",iostat=ierr, &
                       iomsg=error_message)
    call io_error(pname,fname,1,__LINE__,__FILE__,ierr,error_message)
    !
    if (nr == 2) then
      !
      write (inrm,20) nfqpts,nface
      !
      allocate ( ipface(1:2,1:nface) , source=0 , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"ipface",1,__LINE__,__FILE__,ierr,error_message)
      !
    else if (nr == 3) then
      !
      write (inrm,21) nfqpts,nface
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
    n0 = cell(nc)%beg_sp-1
    np = cell(nc)%end_sp-n0
    !
    do k = 1,np
      do l = 1,nr
        txyz(k,l) = xyz(l,n0+k)
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
      qdtr_order = cell(nc)%face(ncf)%qdtr_order
      !
      nqp = geom_solpts(face_geom,qdtr_order)
      !
      ! Interpolate the coordinates of the solution points to get the
      ! coordinates of the flux points on the current face.
      !
      do l = 1,nr
        !
        do kf = 1,nqp
          !
          lf = cell(nc)%face(ncf)%qp_offset + fp%geord(face_geom,qdtr_order)% &
                                            rot(cell(nc)%face(ncf)%rotation)% &
                                              pts(kf)
          !
          val = zero
          do i = interp(left_geom,left_order)%toFace(qdtr_order)%row_ptr(lf)+1,&
                 interp(left_geom,left_order)%toFace(qdtr_order)%row_ptr(lf+1)
            !
            j = interp(left_geom,left_order)%toFace(qdtr_order)%row_idx(i)
            a = interp(left_geom,left_order)%toFace(qdtr_order)%row_val(i)
            !
            val = val + a * txyz(j,l)
            !
          end do
          !
          xyz_qdtr(kf,l) = val
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
      ndim = size(deriv(face_geom,qdtr_order)%d)
      !
      do k = 1,nqp
        !
        do n = 1,ndim
          !
          ! If the current cell face is a x-min, y-max, or z-min face, we need
          ! to reverse the order of the derivatives so that the cross_product
          ! computed below results in a outward facing normal for the left cell
          ! on the face.
          !
          m = n
          if (face_geom /= Geom_Edge) then
            if (any(ncf == [Xmin3D_Face,Ymax3D_Face,Zmin3D_Face])) then
              m = ndim+1 - n
            end if
          end if
          !
          do l = 1,nr
            val = zero
            do i = deriv(face_geom,qdtr_order)%d(n)%row_ptr(k)+1, &
                   deriv(face_geom,qdtr_order)%d(n)%row_ptr(k+1)
              !
              j = deriv(face_geom,qdtr_order)%d(n)%row_idx(i)
              a = deriv(face_geom,qdtr_order)%d(n)%row_val(i)
              !
              val = val + a * xyz_qdtr(j,l)
              !
            end do
            dvec(l,m) = val
          end do
          !
        end do
        !
        ! Finally, the length/area normal to the current face at the current
        ! flux point is just the cross product of the two vectors in dvec
        !
       !dvec(:,3) = cross_product( dvec(:,1) , dvec(:,2) )
        if (face_geom == Geom_Edge) then
          dvec(1,3) = dvec(2,1) - dvec(2,2)
          dvec(2,3) = dvec(1,2) - dvec(1,1)
          dvec(3,3) = zero
        else
          dvec(1,3) = dvec(2,1)*dvec(3,2) - dvec(3,1)*dvec(2,2)
          dvec(2,3) = dvec(3,1)*dvec(1,2) - dvec(1,1)*dvec(3,2)
          dvec(3,3) = dvec(1,1)*dvec(2,2) - dvec(2,1)*dvec(1,2)
        end if
        !
        ! Zero/chop all the normals less than 10*roundoff
        !
        geofq(nf)%v(1:nr,k) = chop( real(dvec(1:nr,3),kind=qp)/qten ) * ten
        !
        if (output_face_normals) then
          !
          npcount = npcount + 1
          !
          if (nr == 2) then
            !
            if (k == 1) then
              ipface(1,nf) = npcount
            else if (k == nqp) then
              ipface(2,nf) = npcount
            end if
            !
            write (inrm,22) (xyz_qdtr(k,l),l=1,nr), &
                            (geofq(nf)%v(l,k),l=1,nr)!,npcount
            !
          else if (nr == 3) then
            !
            if (k == 1) then
              ipface(1,nf) = npcount
            else if (k == qdtr_order+1) then
              ipface(2,nf) = npcount
            else if (k == nqp) then
              ipface(3,nf) = npcount
            else if (face_geom == Geom_Quad) then
              if (k == (qdtr_order+1)*qdtr_order) then
                ipface(4,nf) = npcount
              end if
            end if
            !
            if (face_geom == Geom_Tria) ipface(4,nf) = ipface(3,nf)
            !
            write (inrm,23) (xyz_qdtr(k,l),l=1,nr), &
                            (geofq(nf)%v(l,k),l=1,nr)!,npcount
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
  ! Copy geofq into facenrm_qdtpts
  !
  allocate ( facenrm_qdtpts(1:nr,1:n_flx_qdtpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"facenrm_qdtpts",1,__LINE__,__FILE__, &
                   ierr,error_message)
  !$ACC UPDATE DEVICE(facenrm_qdtpts)
  !
  n = 0
  do nf = 1,size(geofq)
    do k = 1,size(geofq(nf)%v,dim=2)
      n = n + 1
      facenrm_qdtpts(1:nr,n) = geofq(nf)%v(1:nr,k)
    end do
  end do
  !$ACC UPDATE DEVICE(facenrm_qdtpts)
  !
  if (allocated(geofq)) then
    deallocate ( geofq , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"geofq",2,__LINE__,__FILE__,ierr,error_message)
  end if
  !
  !
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
end subroutine compute_quadrature_face_normals
!
!###############################################################################
!
subroutine create_time_metrics()
  !
  !.. Use Statements ..
  use geovar,           only : nr,ncell,cell,n_solpts
  use element_type_mod, only : element_t,std_elem
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,n,nc,ni,n0,np,nd,ierr
  integer  :: this_geom,this_order
  real(wp) :: dx,dy,dz,jacobian
  real(wp) :: d1,d2,d3
  real(wp) :: d1dx,d2dx,d3dx
  real(wp) :: d1dy,d2dy,d3dy
  real(wp) :: d1dz,d2dz,d3dz
  !
  !.. Local Arrays ..
  real(wp) :: del(1:nr)
  !
  !.. Local Allocatable Arrays ..
  type(element_t), allocatable :: del_ref_coord(:,:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_time_metrics"
  character(len=*), parameter :: char_ijk(1:3) = ["xi  ","eta ","zeta"]
  !
  logical(lk), parameter :: account_for_boundaries = true
 !logical(lk), parameter :: account_for_boundaries = fals
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  allocate ( del_ref_coord , source=std_elem , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"del_ref_coord",1,__LINE__,__FILE__,ierr,error_message)
  !
  do this_order = lbound(std_elem,dim=2),ubound(std_elem,dim=2)
    do this_geom = lbound(std_elem,dim=1),ubound(std_elem,dim=1)
      !
      if (.not. allocated(del_ref_coord(this_geom,this_order)%pts)) cycle
      !
      nd = size(std_elem(this_geom,this_order)%pts,dim=1)
      np = size(std_elem(this_geom,this_order)%pts,dim=2)
      !
      do n = 1,np
        !
        do l = 1,nd
          dx = huge(zero)
          do i = 1,np
            if (i == n) cycle
            dy = abs( std_elem(this_geom,this_order)%pts(l,n) - &
                      std_elem(this_geom,this_order)%pts(l,i) )
            if (dy > eps6) dx = min(dx,dy)
          end do
          if (account_for_boundaries) then
            dy = abs( std_elem(this_geom,this_order)%pts(l,n) - one )
            if (dy > eps6) dx = min(dx,dy)
            if (this_geom == Geom_Tria .and. .not. tria_is_biunit) then
              dy = abs( std_elem(this_geom,this_order)%pts(l,n) - zero )
              if (dy > eps6) dx = min(dx,dy)
            else
              dy = abs( std_elem(this_geom,this_order)%pts(l,n) + one )
              if (dy > eps6) dx = min(dx,dy)
            end if
          end if
          del_ref_coord(this_geom,this_order)%pts(l,n) = one/dx
        end do
        !
      end do
      !
      write (140,100) Geom_Name(this_geom),this_order
      do l = 1,nd
        ni = this_order+1
        write (140,102) trim(adjustl(char_ijk(l)))
        if (nd == 1) then
          write (140,101) &
            (one/del_ref_coord(this_geom,this_order)%pts(l,i),i=1,ni)
        else if (nd == 2) then
          do j = ni,1,-1
            write (140,101) &
              (one/del_ref_coord(this_geom,this_order)%pts(l,(j-1)*ni+i),i=1,ni)
          end do
        else if (nd == 3) then
          do k = 1,ni
            do j = ni,1,-1
              write (140,101) &
  (one/del_ref_coord(this_geom,this_order)%pts(l,(k-1)*ni*ni+(j-1)*ni+i),i=1,ni)
            end do
            write (140,*)
          end do
        end if
      end do
      100 format (/,"Geom = ",a," ; Order = ",i0)
      101 format (2x,100es14.6)
      102 format (2x,"delta_",a)
      !
    end do
  end do
  !
  allocate ( time_metrics(1:nr+1,1:n_solpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"time_metrics",1,__LINE__,__FILE__,ierr,error_message)
  !
  do nc = 1,ncell
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    n0 = cell(nc)%beg_sp-1
    np = cell(nc)%end_sp-n0
    !
   !do n = 1,np
   !  !
   !  jacobian = jac(n0+n)
   !  !
   !  del(1:nr) = del_ref_coord(this_geom,this_order)%pts(1:nr,n)
   !  !
   !  do l = 1,nr
   !    !
   !    dx = jacobian * dot( del(1:nr) , metrics_array(n0+n,1:nr,l) )
   !    !
   !    time_metrics(l,n0+n) = one/dx
   !    !
   !  end do
   !  !
   !  time_metrics(nr+1,n0+n) = norm2(time_metrics(1:nr,n0+n))
   !  !
   !end do
    !
    if (nr == 3) then
      !
      do n = 1,np
        !
        d1 = del_ref_coord(this_geom,this_order)%pts(1,n)
        d2 = del_ref_coord(this_geom,this_order)%pts(2,n)
        d3 = del_ref_coord(this_geom,this_order)%pts(3,n)
        !
        jacobian = jac(n0+n)
        !
        d1dx = metrics_array(n0+n,1,1)*jacobian
        d2dx = metrics_array(n0+n,2,1)*jacobian
        d3dx = metrics_array(n0+n,3,1)*jacobian
        !
        d1dy = metrics_array(n0+n,1,2)*jacobian
        d2dy = metrics_array(n0+n,2,2)*jacobian
        d3dy = metrics_array(n0+n,3,2)*jacobian
        !
        d1dz = metrics_array(n0+n,1,3)*jacobian
        d2dz = metrics_array(n0+n,2,3)*jacobian
        d3dz = metrics_array(n0+n,3,3)*jacobian
        !
        dx = one / (d1*d1dx + d2*d2dx + d3*d3dx)
        dy = one / (d1*d1dy + d2*d2dy + d3*d3dy)
        dz = one / (d1*d1dz + d2*d2dz + d3*d3dz)
        !
        time_metrics(1,n0+n) = dx
        time_metrics(2,n0+n) = dy
        time_metrics(3,n0+n) = dz
        time_metrics(4,n0+n) = sqrt(dx*dx + dy*dy + dz*dz)
        !
      end do
      !
    else if (nr == 2) then
      !
      do n = 1,np
        !
        d1 = del_ref_coord(this_geom,this_order)%pts(1,n)
        d2 = del_ref_coord(this_geom,this_order)%pts(2,n)
        !
        jacobian = jac(n0+n)
        !
        d1dx = metrics_array(n0+n,1,1)*jacobian
        d2dx = metrics_array(n0+n,2,1)*jacobian
        !
        d1dy = metrics_array(n0+n,1,2)*jacobian
        d2dy = metrics_array(n0+n,2,2)*jacobian
        !
        dx = one / (d1*d1dx + d2*d2dx)
        dy = one / (d1*d1dy + d2*d2dy)
        !
        time_metrics(1,n0+n) = dx
        time_metrics(2,n0+n) = dy
        time_metrics(3,n0+n) = sqrt(dx*dx + dy*dy)
        !
      end do
      !
    else
      !
      do n = 1,np
        !
        d1 = del_ref_coord(this_geom,this_order)%pts(1,n)
        !
        jacobian = jac(n0+n)
        !
        d1dx = metrics_array(n0+n,1,1)*jacobian
        !
        dx = one / (d1*d1dx)
        !
        time_metrics(1,n0+n) = dx
        time_metrics(2,n0+n) = sqrt(dx*dx)
        !
      end do
      !
    end if
    !
  end do
  !
  deallocate ( del_ref_coord , stat=ierr , errmsg=error_message )
  call alloc_error(pname,"del_ref_coord",2,__LINE__,__FILE__,ierr,error_message)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine create_time_metrics
!
!###############################################################################
!
pure function compute_cell_volume(this,geom,order) result(return_value)
  !
  !.. Use Statements ..
  use element_type_mod, only : std_elem
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
pure function compute_face_area(this,geom,order) result(return_value)
  !
  !.. Use Statements ..
  use element_type_mod, only : std_elem
  !
  !.. Passed Argument from Invoking Object ..
  class(face_nrm_t), intent(in) :: this
  !
  !.. Formal Arguments ..
  integer, intent(in) :: geom
  integer, intent(in) :: order
  !
  !.. Function Result ..
  real(wp) :: return_value
  !
  !.. Local Scalars ..
  integer :: n
  !
continue
  !
  ! Sum the product of the Jacobian and quadrature weight for each point
  ! in this cell to compute the volume of the current cell
  !
  do n = 1,size(this%v,dim=2)
    return_value = return_value + norm2(this%v(:,n)) * &
                                  std_elem(geom,order)%wts(n)
  end do
  !
end function compute_face_area
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
  integer  :: i,j,k,l,m,n,nr,np,ierr
  real(wp) :: a,val
  character(len=50) :: array_name
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "create_cell_metrics"
  !
continue
  !
  np = size(txyz,dim=1)
  nr = size(txyz,dim=2)
  !
  ! Compute the metrics for all points in this cell
  !
  ! USING F2003 AUTO-REALLOCATION
 !this%met = deriv(geom,order)%grad( txyz(:,:) )
  allocate ( this%met(1:np,1:nr,1:nr) , source=zero , &
             stat=ierr , errmsg=error_message )
  write (array_name,1) cell_idx,"met"
  call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
  !
  do l = 1,nr
    do m = 1,nr
      do k = 1,np
        val = zero
        do i = deriv(geom,order)%d(l)%row_ptr(k)+1, &
               deriv(geom,order)%d(l)%row_ptr(k+1)
          !
          j = deriv(geom,order)%d(l)%row_idx(i)
          a = deriv(geom,order)%d(l)%row_val(i)
          !
          val = val + a * txyz(j,m)
          !
        end do
        this%met(k,m,l) = val
      end do
    end do
  end do
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
  write (array_name,1) cell_idx,"jac"
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
    !
   !! USING F2003 AUTO-REALLOCATION
   !this%dJ = deriv(geom,order)%grad( this%jac(:) )
    allocate ( this%dJ(1:np,1:nr) , source=zero , &
               stat=ierr , errmsg=error_message )
    write (array_name,1) cell_idx,"dJ"
    call alloc_error(pname,array_name,1,__LINE__,__FILE__,ierr,error_message)
    !
    do l = 1,nr
      do k = 1,np
        val = zero
        do i = deriv(geom,order)%d(l)%row_ptr(k)+1, &
               deriv(geom,order)%d(l)%row_ptr(k+1)
          !
          j = deriv(geom,order)%d(l)%row_idx(i)
          a = deriv(geom,order)%d(l)%row_val(i)
          !
          val = val + a * this%jac(j)
          !
        end do
        this%dJ(k,l) = val
      end do
    end do
   !do ref_dir = 1,nr
   !  dJ(ref_dir,n1:n2) = deriv(this_geom,this_order)% &
   !                                       d(ref_dir)% &
   !                      mv_mult( jacobian(n1:n2) )
   !end do
  end if
  !
  ! Format Statements
  !
  1 format ("metrics_dt(",i0,")%",a)
  !
end subroutine create_cell_metrics
!
!###############################################################################
!
subroutine compute_cell_metrics_3d(met,jcb,wjac)
  !
  !.. Use Statements ..
  use geovar,            only : ncell,cell,xyz
  use element_type_mod,  only : std_elem
  use derivatives_mod,   only : deriv
  !
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  use omp_lib, only : omp_get_thread_num
  use omp_lib, only : omp_get_num_threads
#endif
#endif
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: met(:,:,:)
  real(wp), intent(inout) :: jcb(:)
  real(wp), intent(inout) :: wjac(:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,jacob
  real(wp) :: dxdi,dydi,dzdi
  real(wp) :: dxdj,dydj,dzdj
  real(wp) :: dxdk,dydk,dzdk
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_cell_metrics_3d"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  !$OMP PARALLEL &
  !$OMP&         DEFAULT(SHARED) &
  !$OMP&         PRIVATE(nc,n0,np,this_geom,this_order,k,i,n,a,jacob,tid, &
  !$OMP&                 dxdi,dydi,dzdi,dxdj,dydj,dzdj,dxdk,dydk,dzdk)
  !
  tid = omp_get_thread_num()
  nthreads = omp_get_num_threads()
  !
  if (mypnum == glb_root) then
    if (tid == glb_root) then
      write (iout,thread_report) pname,nthreads
    end if
  end if
  !
  chunk_size = ceiling( real(ncell,kind=wp) / real(nthreads,kind=wp) )
  chunk_size = chunk_size / 10
  !
  !$OMP DO SCHEDULE(RUNTIME)
 !!$OMP DO SCHEDULE(GUIDED,chunk_size)
 !!$OMP DO SCHEDULE(STATIC,chunk_size)
#endif
#endif
  metrics_loop: do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Compute the metrics and jacobian for this cell
    !
    do k = 1,np
      !
      ! Xi/I derivatives
      !
      dxdi = zero
      dydi = zero
      dzdi = zero
      !
      do i = deriv(this_geom,this_order)%d(1)%row_ptr(k)+1, &
             deriv(this_geom,this_order)%d(1)%row_ptr(k+1)
        !
        n = deriv(this_geom,this_order)%d(1)%row_idx(i)
        a = deriv(this_geom,this_order)%d(1)%row_val(i)
        !
        dxdi = dxdi + a * xyz(1,n0+n)
        dydi = dydi + a * xyz(2,n0+n)
        dzdi = dzdi + a * xyz(3,n0+n)
        !
      end do
      !
      ! Eta/J derivatives
      !
      dxdj = zero
      dydj = zero
      dzdj = zero
      !
      do i = deriv(this_geom,this_order)%d(2)%row_ptr(k)+1, &
             deriv(this_geom,this_order)%d(2)%row_ptr(k+1)
        !
        n = deriv(this_geom,this_order)%d(2)%row_idx(i)
        a = deriv(this_geom,this_order)%d(2)%row_val(i)
        !
        dxdj = dxdj + a * xyz(1,n0+n)
        dydj = dydj + a * xyz(2,n0+n)
        dzdj = dzdj + a * xyz(3,n0+n)
        !
      end do
      !
      ! Zeta/K derivatives
      !
      dxdk = zero
      dydk = zero
      dzdk = zero
      !
      do i = deriv(this_geom,this_order)%d(3)%row_ptr(k)+1, &
             deriv(this_geom,this_order)%d(3)%row_ptr(k+1)
        !
        n = deriv(this_geom,this_order)%d(3)%row_idx(i)
        a = deriv(this_geom,this_order)%d(3)%row_val(i)
        !
        dxdk = dxdk + a * xyz(1,n0+n)
        dydk = dydk + a * xyz(2,n0+n)
        dzdk = dzdk + a * xyz(3,n0+n)
        !
      end do
      !
      jacob = dxdk*(dydi*dzdj - dydj*dzdi) + &
              dydk*(dxdj*dzdi - dxdi*dzdj) + &
              dzdk*(dxdi*dydj - dxdj*dydi)
      !
      met(n0+k,1,1) = ( dydj*dzdk - dydk*dzdj ) / jacob
      met(n0+k,2,1) = ( dydk*dzdi - dydi*dzdk ) / jacob
      met(n0+k,3,1) = ( dydi*dzdj - dydj*dzdi ) / jacob
      !
      met(n0+k,1,2) = ( dxdk*dzdj - dxdj*dzdk ) / jacob
      met(n0+k,2,2) = ( dxdi*dzdk - dxdk*dzdi ) / jacob
      met(n0+k,3,2) = ( dxdj*dzdi - dxdi*dzdj ) / jacob
      !
      met(n0+k,1,3) = ( dxdj*dydk - dxdk*dydj ) / jacob
      met(n0+k,2,3) = ( dxdk*dydi - dxdi*dydk ) / jacob
      met(n0+k,3,3) = ( dxdi*dydj - dxdj*dydi ) / jacob
      !
      jcb(n0+k) = jacob
      !
      wjac(n0+k) = jacob * std_elem(this_geom,this_order)%wts(k)
      !
    end do
    !
  end do metrics_loop
#ifdef DEBUG_ON
#ifdef ENABLE_OPENMP
  !$OMP END DO
  !$OMP END PARALLEL
#endif
#endif
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_cell_metrics_3d
!
!###############################################################################
!
subroutine compute_cell_metrics_2d(met,jcb,wjac)
  !
  !.. Use Statements ..
  use geovar,            only : ncell,cell,xyz
  use element_type_mod,  only : std_elem
  use derivatives_mod,   only : deriv
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: met(:,:,:)
  real(wp), intent(inout) :: jcb(:)
  real(wp), intent(inout) :: wjac(:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,jacob
  real(wp) :: dxdi,dydi
  real(wp) :: dxdj,dydj
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_cell_metrics_2d"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  metrics_loop: do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Compute the metrics and jacobian for this cell
    !
    do k = 1,np
      !
      ! Xi/I derivatives
      !
      dxdi = zero
      dydi = zero
      !
      do i = deriv(this_geom,this_order)%d(1)%row_ptr(k)+1, &
             deriv(this_geom,this_order)%d(1)%row_ptr(k+1)
        !
        n = deriv(this_geom,this_order)%d(1)%row_idx(i)
        a = deriv(this_geom,this_order)%d(1)%row_val(i)
        !
        dxdi = dxdi + a * xyz(1,n0+n)
        dydi = dydi + a * xyz(2,n0+n)
        !
      end do
      !
      ! Eta/J derivatives
      !
      dxdj = zero
      dydj = zero
      !
      do i = deriv(this_geom,this_order)%d(2)%row_ptr(k)+1, &
             deriv(this_geom,this_order)%d(2)%row_ptr(k+1)
        !
        n = deriv(this_geom,this_order)%d(2)%row_idx(i)
        a = deriv(this_geom,this_order)%d(2)%row_val(i)
        !
        dxdj = dxdj + a * xyz(1,n0+n)
        dydj = dydj + a * xyz(2,n0+n)
        !
      end do
      !
      jacob = dxdi*dydj - dxdj*dydi
      !
      met(n0+k,1,1) =  dydj / jacob
      met(n0+k,2,1) = -dydi / jacob
      !
      met(n0+k,1,2) = -dxdj / jacob
      met(n0+k,2,2) =  dxdi / jacob
      !
      jcb(n0+k) = jacob
      !
      wjac(n0+k) = jacob * std_elem(this_geom,this_order)%wts(k)
      !
    end do
    !
  end do metrics_loop
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_cell_metrics_2d
!
!###############################################################################
!
subroutine compute_cell_metrics_1d(met,jcb,wjac)
  !
  !.. Use Statements ..
  use geovar,            only : ncell,cell,xyz
  use element_type_mod,  only : std_elem
  use derivatives_mod,   only : deriv
  !
  !.. Formal Arguments ..
  real(wp), intent(inout) :: met(:,:,:)
  real(wp), intent(inout) :: jcb(:)
  real(wp), intent(inout) :: wjac(:)
  !
  !.. Local Scalars ..
  integer  :: i,k,n,nc,np,n0
  integer  :: this_geom,this_order
  real(wp) :: a,dxdi
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_cell_metrics_1d"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  metrics_loop: do nc = 1,ncell
    !
    n0 = cell(nc)%beg_sp - 1
    np = cell(nc)%end_sp - n0
    !
    this_geom = cell(nc)%geom
    this_order = cell(nc)%order
    !
    ! Compute the metrics and jacobian for this cell
    !
    do k = 1,np
      !
      ! Xi/I derivatives
      !
      dxdi = zero
      !
      do i = deriv(this_geom,this_order)%d(1)%row_ptr(k)+1, &
             deriv(this_geom,this_order)%d(1)%row_ptr(k+1)
        !
        n = deriv(this_geom,this_order)%d(1)%row_idx(i)
        a = deriv(this_geom,this_order)%d(1)%row_val(i)
        !
        dxdi = dxdi + a * xyz(1,n0+n)
        !
      end do
      !
      met(n0+k,1,1) = one / dxdi
      !
      jcb(n0+k) = dxdi
      !
      wjac(n0+k) = dxdi * std_elem(this_geom,this_order)%wts(k)
      !
    end do
    !
  end do metrics_loop
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_cell_metrics_1d
!
!###############################################################################
!
subroutine compute_volume_statistics(min_ncv,max_ncv,vol_min,vol_max, &
                                     min_jcl,max_jcl,min_jsp,max_jsp, &
                                     jac_min,jac_max,jcb,wjac)
  !
  !.. Use Statements ..
  use geovar, only : ncell,cell
  use ovar,   only : local_volume
  !
  !.. Formal Arguments ..
  integer,  intent(inout) :: min_ncv,max_ncv
  integer,  intent(inout) :: min_jcl,max_jcl
  integer,  intent(inout) :: min_jsp,max_jsp
  real(wp), intent(inout) :: vol_min,vol_max
  real(wp), intent(inout) :: jac_min,jac_max
  real(wp),    intent(in) :: jcb(:)
  real(wp),    intent(in) :: wjac(:)
  !
  !.. Local Scalars ..
  integer  :: n,nc,n1,n2
  !
  !.. Local Allocatable Arrays ..
  real(wp), allocatable :: vol(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "compute_volume_statistics"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  !
  ! Compute the volume of each cell
  !
  vol = init_array(ncell,zero) ! F2003 AUTO-REALLOC
  !
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    !
    vol(nc) = sum( wjac(n1:n2) )
    !
  end do
  !
  ! Compute the total volume of the local partition
  !
  local_volume = sum( vol(1:ncell) )
  !
  ! Get the location and value of the minimum and maximum cell volume
  !
  min_ncv = minloc( vol(1:ncell) , dim=1 )
  max_ncv = maxloc( vol(1:ncell) , dim=1 )
  !
  vol_min = vol(min_ncv)
  vol_max = vol(max_ncv)
  !
  ! Get the location and value of the minimum and maximum Jacobian
  !
  min_jcl = 1
  max_jcl = 1
  !
  min_jsp = 1
  max_jsp = 1
  !
  jac_min = huge(zero)
  jac_max = -huge(zero)
  !
  do nc = 1,ncell
    !
    n1 = cell(nc)%beg_sp
    n2 = cell(nc)%end_sp
    !
    do n = n1,n2
      !
      if (jcb(n) < jac_min) then
        jac_min = jcb(n)
        min_jcl = nc
        min_jsp = n-n1+1
      end if
      !
      if (jcb(n) > jac_max) then
        jac_max = jcb(n)
        max_jcl = nc
        max_jsp = n-n1+1
      end if
      !
    end do
    !
  end do
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine compute_volume_statistics
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
