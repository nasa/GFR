module ovar
  !
  use module_kind_types, only : vbp
  !
  implicit none
  !
  private
  !
 !include "kind_defs.h"
  !
  real(vbp), public, save :: gam
  !
  real(vbp), public, save :: cfl
  real(vbp), public, save :: dt_in
  real(vbp), public, save :: dt_global
  !
  integer,   public, save :: itcur
  real(vbp), public, save :: d_t
  real(vbp), public, save :: time
  !
  integer,   public, save :: itrst
  real(vbp), public, save :: tmrst
  real(vbp), public, save :: dtrst
  !
  real(vbp), public, save, dimension(4) :: bc_in
  real(vbp), public, save, dimension(4) :: bc_top
  real(vbp), public, save, dimension(4) :: bc_left
  real(vbp), public, save, dimension(4) :: bc_right
  real(vbp), public, save, dimension(3) :: bc_invar
  !
  real(vbp), public, save :: p_exit
  real(vbp), public, save :: local_volume
  real(vbp), public, save :: total_volume
  !
  real(vbp), public, save :: betark(4,4)
  !
  integer, public, save :: itestcase
  integer, public, save :: iop_solpts
  integer, public, save :: iop_flxpts
  integer, public, save :: iop_tripts
  integer, public, save :: iop_qdtpts
  integer, public, save :: irestart
  integer, public, save :: flux_divergence
  integer, public, save :: iop_correction
  integer, public, save :: num_ts
  integer, public, save :: flux_method
  integer, public, save :: init_loc
  integer, public, save :: error_loc
  integer, public, save :: dflux_opt
  !
  integer,   public, save :: metis_option
  !
  integer,   public, save :: Filter_Option
  integer,   public, save :: Filter_Cutoff
  integer,   public, save :: Filter_Order
  real(vbp), public, save :: Filter_Alpha
  !
  integer,   public, save :: nloc
  integer,   public, save :: kloc
  real(vbp), public, save :: rlinf
  real(vbp), public, save :: rlone
  real(vbp), public, save :: rltwo
  !
  integer,   public, save, allocatable, dimension(:) :: unloc
  integer,   public, save, allocatable, dimension(:) :: ukloc
  real(vbp), public, save, allocatable, dimension(:) :: ulinf
  real(vbp), public, save, allocatable, dimension(:) :: ulone
  real(vbp), public, save, allocatable, dimension(:) :: ultwo
  !
  integer,   public, save, allocatable, dimension(:) :: vnloc
  integer,   public, save, allocatable, dimension(:) :: vkloc
  real(vbp), public, save, allocatable, dimension(:) :: vlinf
  real(vbp), public, save, allocatable, dimension(:) :: vlone
  real(vbp), public, save, allocatable, dimension(:) :: vltwo
  !
  real(vbp), public, save, dimension(1:3,1:3,1:5) :: advec_error
  real(vbp), public, save, dimension(1:5,1:2,1:5) :: vortex_error
  !
  integer,   public, save :: nentmax
  integer,   public, save :: kentmax
  real(vbp), public, save :: xentmax
  real(vbp), public, save :: yentmax
  real(vbp), public, save :: entmax
  real(vbp), public, save :: entL2
  !
  integer, public, save :: nst
  integer, public, save :: num_stages
  integer, public, save :: iRKtype
  !
  integer, public, parameter :: Classic_RK = 1
  integer, public, parameter :: TVD2_RK    = 2
  integer, public, parameter :: TVD3_RK    = 3
  integer, public, parameter :: TVD4_RK    = 4
  integer, public, parameter :: CK4_RK     = 5
  !
  integer, public, save :: output_solution
  integer, public, save :: output_interval
  integer, public, save :: restart_interval
  integer, public, save :: results_interval
  !
  real(vbp), public, save :: rhoref
  real(vbp), public, save :: tref
  real(vbp), public, save :: pref
  real(vbp), public, save :: aref
  real(vbp), public, save :: machref
  real(vbp), public, save :: aoaref
  real(vbp), public, save :: uref
  real(vbp), public, save :: vxref
  real(vbp), public, save :: vyref
  real(vbp), public, save :: rgasref
  !
  logical,   public, save :: bound_error
  real(vbp), public, save :: error_bound_limit
  !
  real(vbp), public, save :: Final_Time
  integer,   public, save :: Total_Periods
  !
  logical,   public, save :: walltime_expiring
  real(vbp), public, save :: max_wall_time
  real(vbp), public, save :: wall_time_left
  real(vbp), public, save :: beg_wall_time
  real(vbp), public, save :: end_wall_time
  real(vbp), public, save :: wall_time
  real(vbp), public, save :: work_units
  !
  ! Using Intel 2013 compilers
  real(vbp), public, parameter :: one_work_unit = 6.9205_vbp
  ! Westmere nodes
 !real(vbp), public, parameter :: one_work_unit = 25.551_vbp/3.0_vbp
  ! Sandy Bridge nodes
 !real(vbp), public, parameter :: one_work_unit = 22.808_vbp/3.0_vbp
  ! Ivy Bridge nodes
 !real(vbp), public, parameter :: one_work_unit = 19.096_vbp/3.0_vbp
  ! Using Intel version 10 compilers
 !real(vbp), public, parameter :: one_work_unit = 7.026_vbp
  !
  real(vbp), public, save :: Period_Time
  integer,   public, save :: Period_Timesteps
  real(vbp), public, save :: Grid_Max_X
  real(vbp), public, save :: Grid_Max_Y
  real(vbp), public, save :: Grid_Min_X
  real(vbp), public, save :: Grid_Min_Y
  real(vbp), public, save :: Grid_Length_X
  real(vbp), public, save :: Grid_Length_Y
  !
  real(vbp), public, save :: Density_Minimum_Residual
  real(vbp), public, save :: Density_Maximum_Residual
  logical,   public, save :: vortex_has_collapsed
  integer,   public, save :: final_countdown
  !
end module ovar
!
!###############################################################################
!
