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
  !
contains
!
!###############################################################################
!
subroutine get_metrics()
  !
  !.. Use Statements ..
  use geovar, only : nr,nface,n_totpts,n_totcel,ncell
  use geovar, only : nodes_of_cell_ptr,nodes_of_cell
  use geovar, only : nodes_on_face,pts_in_cell_ptr
  use geovar, only : xyz_nodes,cell_geom,cell_order
  use geovar, only : cells_on_face
  use geovar, only : solpts_1d,solwts_1d,rs
  !
  use quadrature_mod, only : std_elem
  !
  use ovar, only : local_volume,total_volume
  use ovar, only : using_curved_boundaries
  use ovar, only : use_3rd_derivative_metrics
  !
  use metvar, only : metrics,metrics2,metrics3
  use metvar, only : jacobian,dJ,volj,geofa
  !
  use order_mod, only : n_sp1de,geom_solpts
  !
  !.. Local Scalars ..
  integer  :: i,j,k,l,ierr
  integer  :: n,nc,n1,n2,nf,np,p1,p2
  integer  :: edge_pts,nfnods
  integer  :: this_geom,this_order
  real(wp) :: tvol,tria_coef,tc2,r,s,t
  real(r8) :: loc_vol,tot_vol
  real(wp) :: val,rm1,rp1,sm1,sp1,tm1,tp1,zm1,zp1
  character(len=1000) :: error_message
  !
  !.. Local Arrays ..
  real(wp), dimension(1:nr,1:8) :: nodes
  !
  !.. Local Pointers ..
  real(wp), pointer, contiguous :: weights(:)
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "get_metrics"
  !
continue
  !
  ! Initialize weights to disassociated
  !
  weights => null()
  !
  ! Coefficients for triangle metrics depending on
  ! the definition of the standard triangular element
  !
  tria_coef = one - half*deltafun(tria_is_biunit)
  tc2 = half * tria_coef
  !
  ! Allocate the required metrics related arrays
  !
  allocate ( metrics(1:nr,1:nr,1:n_totpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"metrics",1,ierr,error_message)
  !
  allocate ( jacobian(1:n_totpts) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"jacobian",1,ierr,error_message)
  !
  allocate ( volj(1:n_totcel) , source=zero , &
             stat=ierr , errmsg=error_message )
  call alloc_error(pname,"volj",1,ierr,error_message)
  !
  allocate ( geofa(1:nr,1:n_sp1de,1:nface) , source=zero , stat=ierr )
  call alloc_error(pname,"geofa",1,ierr,error_message)
  !
  ! Allocate some arrays needed for curved boundaries
  !
  if (using_curved_boundaries) then
    !
    allocate ( metrics2(1:nr,1:nr,1:nr,1:n_totpts) , source=zero , &
               stat=ierr , errmsg=error_message )
    call alloc_error(pname,"metrics2",1,ierr,error_message)
    !
    allocate ( dJ(1:nr,1:n_totpts) , stat=ierr , errmsg=error_message )
    call alloc_error(pname,"dJ",1,ierr,error_message)
    !
    if (use_3rd_derivative_metrics) then
      !
      allocate ( metrics3(1:nr,1:nr,1:nr,1:nr,1:n_totpts) , source=zero , &
                 stat=ierr , errmsg=error_message )
      call alloc_error(pname,"metrics3",1,ierr,error_message)
      !
    end if
    !
  end if
  !
  local_volume = zero
  total_volume = zero
  !
  ! Loop over the cells to get the metrics at each solution point
  !
  metrics_loop: do nc = 1,ncell
    !
    n1 = nodes_of_cell_ptr(nc)+1
    n2 = nodes_of_cell_ptr(nc+1)
    !
    p1 = pts_in_cell_ptr(nc)+1
    p2 = pts_in_cell_ptr(nc+1)
    !
    this_geom = cell_geom(nc)
    this_order = cell_order(nc)
    !
    ! Associate weights with the quadrature weights of the
    ! current cells geometry and order
    !
    weights => std_elem(this_geom,this_order)%wts
    !
    ! Load the nodes defining the cell into a temporary array
    !
    do n = n1,n2
      nodes(:,n-n1+1) = xyz_nodes(:,nodes_of_cell(n))
    end do
    !
    ! Initialize the volume for this cell to zero
    !
    tvol = zero
    !
    ! Compute at each solution point (depending on the cell type)
    ! x metrics - dx/dr,dx/ds,d2x/dr2,d2x/ds2,d2x/drds,
    ! y metrics - dy/dr,dy/ds,d2y/dr2,d2y/ds2,d2y/drds,
    ! Jacobians - J,dJ/dr,dJ/ds
    !
    if (cell_geom(nc) == Geom_Edge) then
      !
      ! ###################################
      ! #####  Metrics for edge cell  #####
      ! ###################################
      !
      do n = p1,p2
        !
        np = n-p1+1
        !
        ! First derivative metrics
        !
        metrics(1:nr,1,n) = half * (nodes(1:nr,2) - nodes(1:nr,1))
        !                  _____________________________________
        ! jacobian : J = \/ (dx/dr)**2 + (dy/dr)**2 + (dz/dr)**2
        !
        jacobian(n) = rnorm2( metrics(1:nr,1,n) )
        !
        ! Add the jacobian and quadrature weight to this cells volume
        !
        tvol = tvol + jacobian(n)*solwts_1d(np)
        !
      end do
      !
   else if (cell_geom(nc) == Geom_Tria) then
      !
      ! ########################################
      ! #####  Metrics for triangular cell #####
      ! ########################################
      !
      do n = p1,p2
        !
        np = n-p1+1
        !
        ! First derivative metrics
        !
        metrics(1:nr,1,n) = tria_coef * (nodes(1:nr,2) - nodes(1:nr,1))
        !
        metrics(1:nr,2,n) = tria_coef * (nodes(1:nr,3) - nodes(1:nr,1))
        !
        ! Second derivative metrics
        !
        if (allocated(metrics2)) then
          metrics2(:,:,:,n) = zero
        end if
        !
        ! jacobian : J = ( dx/dr * dy/ds ) - ( dx/ds * dy/dr )
        !
        jacobian(n) = compute_jacobian( metrics(1:nr,1:nr,n) )
        !
        ! First derivative of jacobian
        !
        if (allocated(dJ)) then
          !
          ! dJ/dr = [ ( d2x/dr2 * dy/ds ) + ( dx/dr * d2y/drds ) ] - 
          !         [ ( d2x/drds * dy/dr ) + ( dx/ds * d2y/dr2 ) ]
          !
          dJ(1,n) = metrics2(1,1,1,n)*metrics(2,2,n) + &
                    metrics2(2,1,2,n)*metrics(1,1,n) - &
                    metrics2(1,1,2,n)*metrics(2,1,n) - &
                    metrics2(2,1,1,n)*metrics(1,2,n)
          !
          ! dJ/ds = [ ( d2x/dsdr * dy/ds ) + ( dx/dr * d2y/ds2 ) ] - 
          !         [ ( d2x/ds2 * dy/dr ) + ( dx/ds * d2y/dsdr ) ]
          !
          dJ(2,n) = metrics2(1,1,2,n)*metrics(2,2,n) + &
                    metrics2(2,2,2,n)*metrics(1,1,n) - &
                    metrics2(1,2,2,n)*metrics(2,1,n) - &
                    metrics2(2,1,2,n)*metrics(1,2,n)
          !
        end if
        !
        ! Add the jacobian and quadrature weight to this cells volume
        !
        tvol = tvol + jacobian(n)*weights(np)
        !
      end do
      !
    else if (cell_geom(nc) == Geom_Quad) then
      !
      ! ############################################
      ! #####  Metrics for quadrilateral cell  #####
      ! ############################################
      !
      edge_pts = geom_solpts(Geom_Edge,cell_order(nc))
      !
      do n = p1,p2
        !
        np = n-p1+1
        !
        j = (np-1)/edge_pts + 1
        i = np - (j-1)*edge_pts
        !
        rm1 = solpts_1d(i) - one
        rp1 = solpts_1d(i) + one
        !
        sm1 = solpts_1d(j) - one
        sp1 = solpts_1d(j) + one
        !
        ! First derivative metrics
        !
        metrics(1:nr,1,n) = one4 * ( sm1*(nodes(1:nr,1)-nodes(1:nr,2)) + &
                                     sp1*(nodes(1:nr,3)-nodes(1:nr,4)) )
        !
        metrics(1:nr,2,n) = one4 * ( rm1*(nodes(1:nr,1)-nodes(1:nr,4)) + &
                                     rp1*(nodes(1:nr,3)-nodes(1:nr,2)) )
        !
        ! Second derivative metrics
        !
        if (allocated(metrics2)) then
          !
          metrics2(1:nr,1,2,n) = one4 * ( nodes(1:nr,1) - nodes(1:nr,2) + &
                                          nodes(1:nr,3) - nodes(1:nr,4) )
          !
          metrics2(1:nr,2,1,n) = metrics2(1:nr,1,2,n)
          !
        end if
        !
        ! jacobian : J = ( dx/dr * dy/ds ) - ( dx/ds * dy/dr )
        !
        jacobian(n) = compute_jacobian( metrics(1:nr,1:nr,n) )
        !
        ! First derivative of jacobian
        !
        if (allocated(dJ)) then
          !
          ! dJ/dr = [ ( d2x/dr2 * dy/ds ) + ( dx/dr * d2y/drds ) ] - 
          !         [ ( d2x/drds * dy/dr ) + ( dx/ds * d2y/dr2 ) ]
          !
          dJ(1,n) = metrics2(1,1,1,n)*metrics(2,2,n) + &
                    metrics2(2,1,2,n)*metrics(1,1,n) - &
                    metrics2(1,1,2,n)*metrics(2,1,n) - &
                    metrics2(2,1,1,n)*metrics(1,2,n)
          !
          ! dJ/ds = [ ( d2x/dsdr * dy/ds ) + ( dx/dr * d2y/ds2 ) ] - 
          !         [ ( d2x/ds2 * dy/dr ) + ( dx/ds * d2y/dsdr ) ]
          !
          dJ(2,n) = metrics2(1,1,2,n)*metrics(2,2,n) + &
                    metrics2(2,2,2,n)*metrics(1,1,n) - &
                    metrics2(1,2,2,n)*metrics(2,1,n) - &
                    metrics2(2,1,2,n)*metrics(1,2,n)
          !
        end if
        !
        ! Add the jacobian and quadrature weight to this cells volume
        !
        tvol = tvol + jacobian(n)*weights(np)
        !
      end do
      !
      ! ##########################################################################
      ! ##########################################################################
      ! ##########################################################################
      ! ###############  ||||||||||                   ||||||||||  ################
      ! ###############  VVVVVVVVVV  NEED TO WORK ON  VVVVVVVVVV  ################
      ! ##########################################################################
      ! ##########################################################################
      ! ##########################################################################
      !
    else if (cell_geom(nc) == Geom_Tetr) then
      !
      ! ##########################################
      ! #####  Metrics for tetrahedral cell  #####
      ! ##########################################
      !
      write (error_message,101)
      call stop_recon(stop_mpi,pname,error_message)
      !
    else if (cell_geom(nc) == Geom_Pyra) then
      !
      ! ######################################
      ! #####  Metrics for pyramid cell  #####
      ! ######################################
      !
      write (error_message,102)
      call stop_recon(stop_mpi,pname,error_message)
      !
      ! ##########################################################################
      ! ##########################################################################
      ! ##########################################################################
      ! ###############  AAAAAAAAAA  NEED TO WORK ON  AAAAAAAAAA  ################
      ! ###############  ||||||||||                   ||||||||||  ################
      ! ##########################################################################
      ! ##########################################################################
      ! ##########################################################################
      !
    else if (cell_geom(nc) == Geom_Pris) then
      !
      ! ####################################
      ! #####  Metrics for prism cell  #####
      ! ####################################
      !
      edge_pts = geom_solpts(Geom_Edge,cell_order(nc))
      !
      do n = p1,p2
        !
        np = n-p1+1
        !
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        ! ###############  ||||||||||                   ||||||||||  ################
        ! ###############  VVVVVVVVVV   NEEDS WORK!!!   VVVVVVVVVV  ################
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        !
        k = (np-1)/(edge_pts**2) + 1
        j = (np-1)/edge_pts - (k-1)*edge_pts + 1
        i = np - (j-1)*edge_pts - (k-1)*edge_pts*edge_pts
        !
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        ! ###############  AAAAAAAAAA   NEEDS WORK!!!   AAAAAAAAAA  ################
        ! ###############  ||||||||||                   ||||||||||  ################
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        !
        r = tria_coef * (rs(1,i) + deltafun(tria_is_biunit))
        s = tria_coef * (rs(2,i) + deltafun(tria_is_biunit))
        t = tria_coef * (deltafun(tria_is_biunit) - r - s)
        !
        zm1 = tria_coef * (solpts_1d(j) - one)
        zp1 = tria_coef * (solpts_1d(j) + one)
        !
        metrics(1:nr,1,n) = half * ( zm1*(nodes(1:nr,1) - nodes(1:nr,2)) + &
                                     zp1*(nodes(1:nr,5) - nodes(1:nr,4)) )
        !
        metrics(1:nr,2,n) = half * ( zm1*(nodes(1:nr,1) - nodes(1:nr,3)) + &
                                     zp1*(nodes(1:nr,6) - nodes(1:nr,4)) )
        !
        metrics(1:nr,3,n) = half * ( r*(nodes(1:nr,5) - nodes(1:nr,2)) + &
                                     s*(nodes(1:nr,6) - nodes(1:nr,3)) + &
                                     t*(nodes(1:nr,4) - nodes(1:nr,1)) )
        !
        ! Second derivative metrics
        !
        if (allocated(metrics2)) then
          !
          metrics2(1:nr,2,3,n) = tc2 * ( (nodes(1:nr,1) - nodes(1:nr,3)) + &
                                         (nodes(1:nr,6) - nodes(1:nr,4)) )
          !
          metrics2(1:nr,3,2,n) = tc2 * ( (nodes(1:nr,1) - nodes(1:nr,4)) + &
                                         (nodes(1:nr,6) - nodes(1:nr,3)) )
          !
          metrics2(1:nr,3,1,n) = tc2 * ( (nodes(1:nr,1) - nodes(1:nr,4)) + &
                                         (nodes(1:nr,5) - nodes(1:nr,2)) )
          !
          metrics2(1:nr,1,3,n) = tc2 * ( (nodes(1:nr,1) - nodes(1:nr,2)) + &
                                         (nodes(1:nr,5) - nodes(1:nr,4)) )
          !
        end if
        !
        ! Third derivative metrics
        !
        if (allocated(metrics3)) then
          metrics3(:,:,:,:,n) = zero
        end if
        !
        ! jacobian : J = ( dx/dr * dy/ds ) - ( dx/ds * dy/dr )
        !
        jacobian(n) = compute_jacobian( metrics(1:nr,1:nr,n) )
        !
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        ! ###############  ||||||||||                   ||||||||||  ################
        ! ###############  VVVVVVVVVV  NEED TO WORK ON  VVVVVVVVVV  ################
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        !
        ! First derivative of jacobian
        !
        if (allocated(dJ)) then
          !
          ! dJ/dr = [ ( d2x/dr2 * dy/ds ) + ( dx/dr * d2y/drds ) ] - 
          !         [ ( d2x/drds * dy/dr ) + ( dx/ds * d2y/dr2 ) ]
          !
          dJ(1,n) = metrics2(1,1,1,n)*metrics(2,2,n) + &
                    metrics2(2,1,2,n)*metrics(1,1,n) - &
                    metrics2(1,1,2,n)*metrics(2,1,n) - &
                    metrics2(2,1,1,n)*metrics(1,2,n)
          !
          ! dJ/ds = [ ( d2x/dsdr * dy/ds ) + ( dx/dr * d2y/ds2 ) ] - 
          !         [ ( d2x/ds2 * dy/dr ) + ( dx/ds * d2y/dsdr ) ]
          !
          dJ(2,n) = metrics2(1,1,2,n)*metrics(2,2,n) + &
                    metrics2(2,2,2,n)*metrics(1,1,n) - &
                    metrics2(1,2,2,n)*metrics(2,1,n) - &
                    metrics2(2,1,2,n)*metrics(1,2,n)
          !
        end if
        !
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        ! ###############  AAAAAAAAAA  NEED TO WORK ON  AAAAAAAAAA  ################
        ! ###############  ||||||||||                   ||||||||||  ################
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        !
        ! Add the jacobian and quadrature weight to this cells volume
        !
        tvol = tvol + jacobian(n)*weights(np)
        !
      end do
      !
      !
      write (error_message,103)
      call stop_recon(stop_mpi,pname,error_message)
      !
    else if (cell_geom(nc) == Geom_Hexa) then
      !
      ! #########################################
      ! #####  Metrics for hexahedral cell  #####
      ! #########################################
      !
      edge_pts = geom_solpts(Geom_Edge,cell_order(nc))
      !
      do n = p1,p2
        !
        np = n-p1+1
        !
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        ! ###############  ||||||||||                   ||||||||||  ################
        ! ###############  VVVVVVVVVV DOUBLE CHECK THIS VVVVVVVVVV  ################
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        !
        k = (np-1)/(edge_pts**2) + 1
        j = (np-1)/edge_pts - (k-1)*edge_pts + 1
        i = np - (j-1)*edge_pts - (k-1)*edge_pts*edge_pts
        !
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        ! ###############  AAAAAAAAAA DOUBLE CHECK THIS AAAAAAAAAA  ################
        ! ###############  ||||||||||                   ||||||||||  ################
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        !
        rm1 = solpts_1d(i) - one
        rp1 = solpts_1d(i) + one
        !
        sm1 = solpts_1d(j) - one
        sp1 = solpts_1d(j) + one
        !
        tm1 = solpts_1d(k) - one
        tp1 = solpts_1d(k) + one
        !
        metrics(1:nr,1,n) = one8 * ( sm1*tm1*(nodes(1:nr,2) - nodes(1:nr,1)) + &
                                     sm1*tp1*(nodes(1:nr,5) - nodes(1:nr,6)) + &
                                     sp1*tm1*(nodes(1:nr,4) - nodes(1:nr,3)) + &
                                     sp1*tp1*(nodes(1:nr,7) - nodes(1:nr,8)) )
        !
        metrics(1:nr,2,n) = one8 * ( rm1*tm1*(nodes(1:nr,4) - nodes(1:nr,1)) + &
                                     rm1*tp1*(nodes(1:nr,5) - nodes(1:nr,8)) + &
                                     rp1*tm1*(nodes(1:nr,2) - nodes(1:nr,3)) + &
                                     rp1*tp1*(nodes(1:nr,7) - nodes(1:nr,6)) )
        !
        metrics(1:nr,3,n) = one8 * ( rm1*sm1*(nodes(1:nr,5) - nodes(1:nr,1)) + &
                                     rm1*sp1*(nodes(1:nr,4) - nodes(1:nr,8)) + &
                                     rp1*sm1*(nodes(1:nr,2) - nodes(1:nr,6)) + &
                                     rp1*sp1*(nodes(1:nr,7) - nodes(1:nr,3)) )
        !
        ! Second derivative metrics
        !
        if (allocated(metrics2)) then
          !
          metrics2(1:nr,2,1,n) = one8 * ( tm1*(nodes(1:nr,2) - nodes(1:nr,1)) + &
                                          tm1*(nodes(1:nr,4) - nodes(1:nr,3)) + &
                                          tp1*(nodes(1:nr,5) - nodes(1:nr,6)) + &
                                          tp1*(nodes(1:nr,7) - nodes(1:nr,8)) )
          metrics2(1:nr,1,2,n) = metrics2(1:nr,2,1,n)
          !
          metrics2(1:nr,3,2,n) = one8 * ( rm1*(nodes(1:nr,4) - nodes(1:nr,1)) + &
                                          rm1*(nodes(1:nr,5) - nodes(1:nr,8)) + &
                                          rp1*(nodes(1:nr,2) - nodes(1:nr,3)) + &
                                          rp1*(nodes(1:nr,7) - nodes(1:nr,6)) )
          metrics2(1:nr,2,3,n) = metrics2(1:nr,3,2,n)
          !
          metrics2(1:nr,1,3,n) = one8 * ( sm1*(nodes(1:nr,5) - nodes(1:nr,1)) + &
                                          sm1*(nodes(1:nr,2) - nodes(1:nr,6)) + &
                                          sp1*(nodes(1:nr,4) - nodes(1:nr,8)) + &
                                          sp1*(nodes(1:nr,7) - nodes(1:nr,3)) )
          metrics2(1:nr,3,1,n) = metrics2(1:nr,1,3,n)
          !
        end if
        !
        ! Third derivative metrics
        !
        if (allocated(metrics3)) then
          !
          metrics3(1:nr,3,2,1,n) = one8 * ( nodes(1:nr,2) - nodes(1:nr,1) + &
                                            nodes(1:nr,4) - nodes(1:nr,3) + &
                                            nodes(1:nr,5) - nodes(1:nr,6) + &
                                            nodes(1:nr,7) - nodes(1:nr,8) )
          metrics3(1:nr,2,3,1,n) = metrics3(1:nr,3,2,1,n)
          metrics3(1:nr,3,1,2,n) = metrics3(1:nr,3,2,1,n)
          metrics3(1:nr,1,3,2,n) = metrics3(1:nr,3,2,1,n)
          metrics3(1:nr,2,1,3,n) = metrics3(1:nr,3,2,1,n)
          metrics3(1:nr,1,2,3,n) = metrics3(1:nr,3,2,1,n)
          !
        end if
        !
        ! jacobian : J = ( dx/dr * dy/ds ) - ( dx/ds * dy/dr )
        !
        jacobian(n) = compute_jacobian( metrics(1:nr,1:nr,n) )
        !
        ! First derivative of jacobian
        !
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        ! ###############  ||||||||||                   ||||||||||  ################
        ! ###############  VVVVVVVVVV  NEED TO WORK ON  VVVVVVVVVV  ################
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        !
        if (allocated(dJ)) then
          !
          ! dJ/dr = [ ( d2x/dr2 * dy/ds ) + ( dx/dr * d2y/drds ) ] - 
          !         [ ( d2x/drds * dy/dr ) + ( dx/ds * d2y/dr2 ) ]
          !
          dJ(1,n) = metrics2(1,1,1,n)*metrics(2,2,n) + &
                    metrics2(2,1,2,n)*metrics(1,1,n) - &
                    metrics2(1,1,2,n)*metrics(2,1,n) - &
                    metrics2(2,1,1,n)*metrics(1,2,n)
          !
          ! dJ/ds = [ ( d2x/dsdr * dy/ds ) + ( dx/dr * d2y/ds2 ) ] - 
          !         [ ( d2x/ds2 * dy/dr ) + ( dx/ds * d2y/dsdr ) ]
          !
          dJ(2,n) = metrics2(1,1,2,n)*metrics(2,2,n) + &
                    metrics2(2,2,2,n)*metrics(1,1,n) - &
                    metrics2(1,2,2,n)*metrics(2,1,n) - &
                    metrics2(2,1,2,n)*metrics(1,2,n)
          !
        end if
        !
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        ! ###############  AAAAAAAAAA  NEED TO WORK ON  AAAAAAAAAA  ################
        ! ###############  ||||||||||                   ||||||||||  ################
        ! ##########################################################################
        ! ##########################################################################
        ! ##########################################################################
        !
        ! Add the jacobian and quadrature weight to this cells volume
        !
        tvol = tvol + jacobian(n)*weights(np)
        !
      end do
      !
    else
      !
      ! Unknown or unsupported geometry
      !
      write (error_message,100)
      call stop_recon(stop_mpi,pname,error_message)
      !
    end if
    !
    ! Store the computed volume of the current cell
    !
    volj(nc) = tvol
    !
    ! Add the volume of the current cell to the total volume
    ! of the current local partition
    !
    if (nc <= ncell) local_volume = local_volume + volj(nc)
    !
    ! Disassociate weights before moving to the next cell
    !
    if (associated(weights)) weights => null()
    !
  end do metrics_loop
 !write (*,*) mypnum,local_volume,sum(volj(1:ncell))
  !
  ! ##########################################################################
  ! ##########################################################################
  ! ##########################################################################
  ! ###############  ||||||||||                   ||||||||||  ################
  ! ###############  VVVVVVVVVV  NEED TO WORK ON  VVVVVVVVVV  ################
  ! ##########################################################################
  ! ##########################################################################
  ! ##########################################################################
  !
  ! Get the face normal at the solution points on each face
  !
  ! NOTE: 1. This is correct for the reference triangle defined by
  !          the three nodes [(0,0)-(1,0)-(0,1)] and for the reference
  !          quad defined by the four nodes [(0,0)-(1,0)-(1,1)-(0,1)].
  !       2. This is incorrect (it needs to be multiplied by 0.5) for
  !          the reference triangle defined by the three nodes
  !          [(-1,-1)-(1,-1)-(-1,1)] and for the reference quad defined
  !          by the four nodes [(-1,-1)-(1,-1)-(1,1)-(-1,1)].
  !       3. Since the reference triangle is defined by (1) and the
  !          reference quad is defined by (2), the variable "cl" defined
  !          in the header of the module flux_mod is used to adjust the
  !          face normals between triangles and quads.
  !
  do nf = 1,nface
    !
    ! # of nodes defining this face
    !
    nfnods = geom_nodes(cells_on_face(7,nf))
    !
    ! Coordinates of the nodes on this face
    !
    nodes(1:nr,1:nfnods) = xyz_nodes(1:nr,nodes_on_face(1:nfnods,nf))
    !
    ! !!!!!!!!!!!!!! WARNING: POSSIBLY BAD ASSUMPTIONS MADE HERE !!!!!!!!!!!!!!
    ! Assumption 1: the normal is the same for all flux points on the face
    ! Assumption 2: the cell/face order is uniform throughout the domain
    !               meaning the number of flux points is the same for all
    !               similar face geometries.
    ! !!!!!!!!!!!!!! WARNING: POSSIBLY BAD ASSUMPTIONS MADE HERE !!!!!!!!!!!!!!
    !
    ! Compute the face normal at the first flux point
    !
    geofa(1:nr,1,nf) = face_normal( nodes(1:nr,1:nfnods) )
    !
    ! Copy the face normal from the first flux point to all the other
    ! face flux points (assuming the normal is the same for all of them)
    !
    do n = 2,size(geofa,dim=2)
      geofa(1:nr,n,nf) = geofa(1:nr,1,nf)
    end do
    !
  end do
  !
  ! ##########################################################################
  ! ##########################################################################
  ! ##########################################################################
  ! ###############  AAAAAAAAAA  NEED TO WORK ON  AAAAAAAAAA  ################
  ! ###############  ||||||||||                   ||||||||||  ################
  ! ##########################################################################
  ! ##########################################################################
  ! ##########################################################################
  !
  ! Reduce the volumes of each partition to the root
  ! processor to determine the volume of the global grid.
  !
  tot_vol = real(local_volume,r8)
  loc_vol = real(local_volume,r8)
  !
  if (ncpu > 1) then
    call mpi_allreduce(loc_vol,tot_vol,1,mpi_flttyp, &
                       mpi_sum,mpi_comm_world,mpierr)
  end if
  !
  total_volume = real(tot_vol,kind=wp)
  local_volume = real(loc_vol,kind=wp)
  !
  call mpi_barrier(mpi_comm_world,mpierr)
  !
  ! Format Statements
  !
  100 format (" Encountered an unknown or unsupported cell geometry while", &
              " trying to compute the grid metrics!")
  101 format (" Encountered a tetrahedral grid cell while trying to compute", &
              " the grid metrics! This grid geometry is not currently supported!")
  102 format (" Encountered a pyramid grid cell while trying to compute", &
              " the grid metrics! This grid geometry is not currently supported!")
  103 format (" Encountered a prism grid cell while trying to compute", &
              " the grid metrics! This grid geometry is not currently supported!")
  104 format (" Encountered an invalid metrics derived type!")
  !
end subroutine get_metrics
!
!###############################################################################
!
end module metrics_mod
