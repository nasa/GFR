module ovar
  !
  use module_kind_types, only : wp
  !
  implicit none
  !
  private
  !
 !include "kind_defs.h"
  !
  real(wp), public, save :: gam
  !
  real(wp), public, save :: cfl
  real(wp), public, save :: dt_in
  real(wp), public, save :: dt_global
  !
  integer,   public, save :: itcur
  real(wp), public, save :: d_t
  real(wp), public, save :: time
  !
  integer,   public, save :: itrst
  real(wp), public, save :: tmrst
  real(wp), public, save :: dtrst
  !
  real(wp), public, save, dimension(4) :: bc_in
  real(wp), public, save, dimension(4) :: bc_top
  real(wp), public, save, dimension(4) :: bc_left
  real(wp), public, save, dimension(4) :: bc_right
  real(wp), public, save, dimension(3) :: bc_invar
  !
  real(wp), public, save :: p_exit
  real(wp), public, save :: local_volume
  real(wp), public, save :: total_volume
  !
  real(wp), public, save :: betark(4,4)
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
  real(wp), public, save :: Filter_Alpha
  real(wp), public, save :: Filter_etac
  !
  integer,   public, save :: nloc
  integer,   public, save :: kloc
  real(wp), public, save :: rlinf
  real(wp), public, save :: rlone
  real(wp), public, save :: rltwo
  !
  integer,   public, save, allocatable, dimension(:) :: unloc
  integer,   public, save, allocatable, dimension(:) :: ukloc
  real(wp), public, save, allocatable, dimension(:) :: ulinf
  real(wp), public, save, allocatable, dimension(:) :: ulone
  real(wp), public, save, allocatable, dimension(:) :: ultwo
  !
  integer,   public, save, allocatable, dimension(:) :: vnloc
  integer,   public, save, allocatable, dimension(:) :: vkloc
  real(wp), public, save, allocatable, dimension(:) :: vlinf
  real(wp), public, save, allocatable, dimension(:) :: vlone
  real(wp), public, save, allocatable, dimension(:) :: vltwo
  !
  real(wp), public, save, dimension(1:3,1:3,1:5) :: advec_error
  real(wp), public, save, dimension(1:5,1:2,1:5) :: vortex_error
  !
  integer,   public, save :: nentmax
  integer,   public, save :: kentmax
  real(wp), public, save :: xentmax
  real(wp), public, save :: yentmax
  real(wp), public, save :: entmax
  real(wp), public, save :: entL2
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
  integer, public, save :: iter_out_interval
  !
  real(wp), public, save :: rhoref
  real(wp), public, save :: tref
  real(wp), public, save :: pref
  real(wp), public, save :: aref
  real(wp), public, save :: machref
  real(wp), public, save :: aoaref
  real(wp), public, save :: uref
  real(wp), public, save :: vxref
  real(wp), public, save :: vyref
  real(wp), public, save :: rgasref
  !
  logical,   public, save :: bound_error
  real(wp), public, save :: error_bound_limit
  !
  real(wp), public, save :: Final_Time
  integer,   public, save :: Total_Periods
  !
  logical,   public, save :: walltime_expiring
  real(wp), public, save :: max_wall_time
  real(wp), public, save :: wall_time_left
  real(wp), public, save :: beg_wall_time
  real(wp), public, save :: end_wall_time
  real(wp), public, save :: wall_time
  real(wp), public, save :: work_units
  !
  ! Using Intel 2013 compilers
  real(wp), public, parameter :: one_work_unit = 6.9205_wp
  ! Westmere nodes
 !real(wp), public, parameter :: one_work_unit = 25.551_wp/3.0_wp
  ! Sandy Bridge nodes
 !real(wp), public, parameter :: one_work_unit = 22.808_wp/3.0_wp
  ! Ivy Bridge nodes
 !real(wp), public, parameter :: one_work_unit = 19.096_wp/3.0_wp
  ! Using Intel version 10 compilers
 !real(wp), public, parameter :: one_work_unit = 7.026_wp
  !
  real(wp), public, save :: Period_Time
  integer,   public, save :: Period_Timesteps
  real(wp), public, save :: Grid_Max_X
  real(wp), public, save :: Grid_Max_Y
  real(wp), public, save :: Grid_Min_X
  real(wp), public, save :: Grid_Min_Y
  real(wp), public, save :: Grid_Length_X
  real(wp), public, save :: Grid_Length_Y
  !
  real(wp), public, save :: Density_Minimum_Residual
  real(wp), public, save :: Density_Maximum_Residual
  logical,   public, save :: vortex_has_collapsed
  integer,   public, save :: final_countdown
  !
end module ovar
!
!###############################################################################
!
module ovar_test
  !
  use module_kind_types, only : wp
  use module_kind_types, only : iop_solpts => iop_solpts_ktmod
  use module_kind_types, only : iop_flxpts => iop_flxpts_ktmod
  use module_kind_types, only : iop_qdtpts => iop_qdtpts_ktmod
  use module_kind_types, only : iop_tripts => iop_tripts_ktmod
  use module_kind_types, only : metis_option => metis_opt_ktmod
  use module_kind_types, only : flux_divergence => flux_div_ktmod
  use module_kind_types, only : iop_correction => iop_corr_ktmod
  use module_kind_types, only : dflux_opt => dflux_opt_ktmod
  use module_kind_types, only : Filter_Option => Filter_Option_ktmod
  use module_kind_types, only : Filter_Order => Filter_Order_ktmod
  use module_kind_types, only : Filter_Alpha => Filter_Alpha_ktmod
  use module_kind_types, only : Filter_etac => Filter_etac_ktmod
  use module_kind_types, only : Filter_Cutoff => Filter_Cutoff_ktmod
  use module_kind_types, only : error_bound_limit => err_bnd_lim_ktmod
  use module_kind_types, only : bound_error => bound_error_ktmod
  !
  implicit none
  !
  private
  !
 !include "kind_defs.h"
  !
  public :: iop_solpts
  public :: iop_flxpts
  public :: iop_tripts
  public :: iop_qdtpts
  public :: metis_option
  public :: Filter_Option
  public :: Filter_Cutoff
  public :: Filter_Order
  public :: Filter_Alpha
  public :: Filter_etac
  public :: dflux_opt
  public :: flux_divergence
  public :: iop_correction
  public :: bound_error
  public :: error_bound_limit
  !
  real(wp), public, save :: gam
  !
  real(wp), public, save :: cfl
  real(wp), public, save :: dt_in
  real(wp), public, save :: dt_global
  !
  integer,   public, save :: itcur
  real(wp), public, save :: d_t
  real(wp), public, save :: time
  !
  integer,   public, save :: itrst
  real(wp), public, save :: tmrst
  real(wp), public, save :: dtrst
  !
  real(wp), public, save, dimension(4) :: bc_in
  real(wp), public, save, dimension(4) :: bc_top
  real(wp), public, save, dimension(4) :: bc_left
  real(wp), public, save, dimension(4) :: bc_right
  real(wp), public, save, dimension(3) :: bc_invar
  !
  real(wp), public, save :: p_exit
  real(wp), public, save :: local_volume
  real(wp), public, save :: total_volume
  !
  real(wp), public, save :: betark(4,4)
  !
  integer, public, save :: itestcase
 !integer, public, save :: iop_solpts
 !integer, public, save :: iop_flxpts
 !integer, public, save :: iop_tripts
 !integer, public, save :: iop_qdtpts
  integer, public, save :: irestart
 !integer, public, save :: flux_divergence
 !integer, public, save :: iop_correction
  integer, public, save :: num_ts
  integer, public, save :: flux_method
  integer, public, save :: init_loc
  integer, public, save :: error_loc
 !integer, public, save :: dflux_opt
  !
 !integer,   public, save :: metis_option
  !
 !integer,   public, save :: Filter_Option
 !integer,   public, save :: Filter_Cutoff
 !integer,   public, save :: Filter_Order
 !real(wp), public, save :: Filter_Alpha
 !real(wp), public, save :: Filter_etac
  !
  integer,   public, save :: nloc
  integer,   public, save :: kloc
  real(wp), public, save :: rlinf
  real(wp), public, save :: rlone
  real(wp), public, save :: rltwo
  !
  integer,   public, save, allocatable, dimension(:) :: unloc
  integer,   public, save, allocatable, dimension(:) :: ukloc
  real(wp), public, save, allocatable, dimension(:) :: ulinf
  real(wp), public, save, allocatable, dimension(:) :: ulone
  real(wp), public, save, allocatable, dimension(:) :: ultwo
  !
  integer,   public, save, allocatable, dimension(:) :: vnloc
  integer,   public, save, allocatable, dimension(:) :: vkloc
  real(wp), public, save, allocatable, dimension(:) :: vlinf
  real(wp), public, save, allocatable, dimension(:) :: vlone
  real(wp), public, save, allocatable, dimension(:) :: vltwo
  !
  real(wp), public, save, dimension(1:3,1:3,1:5) :: advec_error
  real(wp), public, save, dimension(1:5,1:2,1:5) :: vortex_error
  !
  integer,   public, save :: nentmax
  integer,   public, save :: kentmax
  real(wp), public, save :: xentmax
  real(wp), public, save :: yentmax
  real(wp), public, save :: entmax
  real(wp), public, save :: entL2
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
  integer, public, save :: iter_out_interval
  !
  real(wp), public, save :: rhoref
  real(wp), public, save :: tref
  real(wp), public, save :: pref
  real(wp), public, save :: aref
  real(wp), public, save :: machref
  real(wp), public, save :: aoaref
  real(wp), public, save :: uref
  real(wp), public, save :: vxref
  real(wp), public, save :: vyref
  real(wp), public, save :: rgasref
  !
 !logical,   public, save :: bound_error
 !real(wp), public, save :: error_bound_limit
  !
  real(wp), public, save :: Final_Time
  integer,   public, save :: Total_Periods
  !
  logical,   public, save :: walltime_expiring
  real(wp), public, save :: max_wall_time
  real(wp), public, save :: wall_time_left
  real(wp), public, save :: beg_wall_time
  real(wp), public, save :: end_wall_time
  real(wp), public, save :: wall_time
  real(wp), public, save :: work_units
  !
  ! Using Intel 2013 compilers
  real(wp), public, parameter :: one_work_unit = 6.9205_wp
  ! Westmere nodes
 !real(wp), public, parameter :: one_work_unit = 25.551_wp/3.0_wp
  ! Sandy Bridge nodes
 !real(wp), public, parameter :: one_work_unit = 22.808_wp/3.0_wp
  ! Ivy Bridge nodes
 !real(wp), public, parameter :: one_work_unit = 19.096_wp/3.0_wp
  ! Using Intel version 10 compilers
 !real(wp), public, parameter :: one_work_unit = 7.026_wp
  !
  real(wp), public, save :: Period_Time
  integer,   public, save :: Period_Timesteps
  real(wp), public, save :: Grid_Max_X
  real(wp), public, save :: Grid_Max_Y
  real(wp), public, save :: Grid_Min_X
  real(wp), public, save :: Grid_Min_Y
  real(wp), public, save :: Grid_Length_X
  real(wp), public, save :: Grid_Length_Y
  !
  real(wp), public, save :: Density_Minimum_Residual
  real(wp), public, save :: Density_Maximum_Residual
  logical,   public, save :: vortex_has_collapsed
  integer,   public, save :: final_countdown
  !
end module ovar_test
!
!###############################################################################
!
