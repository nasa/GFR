module erk_mod
  !
!#define PROFILE_ON
  !
  use module_kind_types, global_debug_timer =>      debug_timer, &
                    pure_global_debug_timer => pure_debug_timer
  use eqn_idx, only : nec,nmx,nmy,nmz,nmb,nme,nee,nq
  use eqn_idx, only : ncs
  !
  use ovar, only : force_nd_subroutines
  !
  implicit none
  !
  private
  !
  public :: update_Classic_RK
  public :: update_CK_RK4
  public :: update_TVD_RK2
  public :: update_TVD_RK3
  public :: update_TVD_RK4
  !
contains
!
!###############################################################################
!
subroutine update_CK_RK4(nst,Apply_Filter,Apply_Limiter)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use flowvar, only : rk_sum,residual,usp,dtsp
  use rhs_mod, only : stabilize_usp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_CK_RK4"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nst == 1) then
    call update_CK_RK4_S1(n_solpts,dtsp,residual,usp,rk_sum)
  else if (nst == 2) then
    call update_CK_RK4_S2(n_solpts,dtsp,residual,usp,rk_sum)
  else if (nst == 3) then
    call update_CK_RK4_S3(n_solpts,dtsp,residual,usp,rk_sum)
  else if (nst == 4) then
    call update_CK_RK4_S4(n_solpts,dtsp,residual,usp,rk_sum)
  else if (nst == 5) then
    call update_CK_RK4_S5(n_solpts,dtsp,residual,usp,rk_sum)
  end if
  !
  ! Apply filtering or limiting to stabilize the solution, and then update
  ! the thermodynamic and transport properties at the solution points
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Limiter,caller_name=pname)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine update_CK_RK4
!
!###############################################################################
!
pure subroutine update_CK_RK4_S1(np,dtsp,residual,usp,rk_sum)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rk_sum(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  real(wp) :: sum_ec,sum_mx,sum_my,sum_mz,sum_ee,sum_eq
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_CK_RK4_S1"
  !
  real(wp), parameter :: b = 1432997174477.0_wp / 9575080441755.0_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = dtsp(n)*residual(n,nec)
      sum_mx = dtsp(n)*residual(n,nmx)
      sum_my = dtsp(n)*residual(n,nmy)
      sum_mz = dtsp(n)*residual(n,nmz)
      sum_ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nmy) = usp(n,nmy) + b*sum_my
      usp(n,nmz) = usp(n,nmz) + b*sum_mz
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nmy) = sum_my
      rk_sum(n,nmz) = sum_mz
      rk_sum(n,nee) = sum_ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = dtsp(n)*residual(n,nec)
      sum_mx = dtsp(n)*residual(n,nmx)
      sum_my = dtsp(n)*residual(n,nmy)
      sum_ee = dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nmy) = usp(n,nmy) + b*sum_my
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nmy) = sum_my
      rk_sum(n,nee) = sum_ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,rk_sum,dtsp,residual)
    do m = 1,nq
      do n = 1,np
        !
        sum_eq = dtsp(n)*residual(n,m)
        !
        usp(n,m) = usp(n,m) + b*sum_eq
        !
        rk_sum(n,m) = sum_eq
        !
      end do
    end do
  end if
  !
end subroutine update_CK_RK4_S1
!
!###############################################################################
!
pure subroutine update_CK_RK4_S2(np,dtsp,residual,usp,rk_sum)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rk_sum(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  real(wp) :: sum_ec,sum_mx,sum_my,sum_mz,sum_ee,sum_eq
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_CK_RK4_S2"
  !
  real(wp), parameter :: a = -567301805773.0_wp / 1357537059087.0_wp
  real(wp), parameter :: b = 5161836677717.0_wp / 13612068292357.0_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = a*rk_sum(n,nec) + dtsp(n)*residual(n,nec)
      sum_mx = a*rk_sum(n,nmx) + dtsp(n)*residual(n,nmx)
      sum_my = a*rk_sum(n,nmy) + dtsp(n)*residual(n,nmy)
      sum_mz = a*rk_sum(n,nmz) + dtsp(n)*residual(n,nmz)
      sum_ee = a*rk_sum(n,nee) + dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nmy) = usp(n,nmy) + b*sum_my
      usp(n,nmz) = usp(n,nmz) + b*sum_mz
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nmy) = sum_my
      rk_sum(n,nmz) = sum_mz
      rk_sum(n,nee) = sum_ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = a*rk_sum(n,nec) + dtsp(n)*residual(n,nec)
      sum_mx = a*rk_sum(n,nmx) + dtsp(n)*residual(n,nmx)
      sum_my = a*rk_sum(n,nmy) + dtsp(n)*residual(n,nmy)
      sum_ee = a*rk_sum(n,nee) + dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nmy) = usp(n,nmy) + b*sum_my
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nmy) = sum_my
      rk_sum(n,nee) = sum_ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,rk_sum,dtsp,residual)
    do m = 1,nq
      do n = 1,np
        !
        sum_eq = a*rk_sum(n,m) + dtsp(n)*residual(n,m)
        !
        usp(n,m) = usp(n,m) + b*sum_eq
        !
        rk_sum(n,m) = sum_eq
        !
      end do
    end do
  end if
  !
end subroutine update_CK_RK4_S2
!
!###############################################################################
!
pure subroutine update_CK_RK4_S3(np,dtsp,residual,usp,rk_sum)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rk_sum(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  real(wp) :: sum_ec,sum_mx,sum_my,sum_mz,sum_ee,sum_eq
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_CK_RK4_S3"
  !
  real(wp), parameter :: a = -2404267990393.0_wp / 2016746695238.0_wp
  real(wp), parameter :: b = 1720146321549.0_wp / 2090206949498.0_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = a*rk_sum(n,nec) + dtsp(n)*residual(n,nec)
      sum_mx = a*rk_sum(n,nmx) + dtsp(n)*residual(n,nmx)
      sum_my = a*rk_sum(n,nmy) + dtsp(n)*residual(n,nmy)
      sum_mz = a*rk_sum(n,nmz) + dtsp(n)*residual(n,nmz)
      sum_ee = a*rk_sum(n,nee) + dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nmy) = usp(n,nmy) + b*sum_my
      usp(n,nmz) = usp(n,nmz) + b*sum_mz
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nmy) = sum_my
      rk_sum(n,nmz) = sum_mz
      rk_sum(n,nee) = sum_ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = a*rk_sum(n,nec) + dtsp(n)*residual(n,nec)
      sum_mx = a*rk_sum(n,nmx) + dtsp(n)*residual(n,nmx)
      sum_my = a*rk_sum(n,nmy) + dtsp(n)*residual(n,nmy)
      sum_ee = a*rk_sum(n,nee) + dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nmy) = usp(n,nmy) + b*sum_my
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nmy) = sum_my
      rk_sum(n,nee) = sum_ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,rk_sum,dtsp,residual)
    do m = 1,nq
      do n = 1,np
        !
        sum_eq = a*rk_sum(n,m) + dtsp(n)*residual(n,m)
        !
        usp(n,m) = usp(n,m) + b*sum_eq
        !
        rk_sum(n,m) = sum_eq
        !
      end do
    end do
  end if
  !
end subroutine update_CK_RK4_S3
!
!###############################################################################
!
pure subroutine update_CK_RK4_S4(np,dtsp,residual,usp,rk_sum)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rk_sum(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  real(wp) :: sum_ec,sum_mx,sum_my,sum_mz,sum_ee,sum_eq
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_CK_RK4_S4"
  !
  real(wp), parameter :: a = -3550918686646.0_wp / 2091501179385.0_wp
  real(wp), parameter :: b = 3134564353537.0_wp / 4481467310338.0_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = a*rk_sum(n,nec) + dtsp(n)*residual(n,nec)
      sum_mx = a*rk_sum(n,nmx) + dtsp(n)*residual(n,nmx)
      sum_my = a*rk_sum(n,nmy) + dtsp(n)*residual(n,nmy)
      sum_mz = a*rk_sum(n,nmz) + dtsp(n)*residual(n,nmz)
      sum_ee = a*rk_sum(n,nee) + dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nmy) = usp(n,nmy) + b*sum_my
      usp(n,nmz) = usp(n,nmz) + b*sum_mz
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nmy) = sum_my
      rk_sum(n,nmz) = sum_mz
      rk_sum(n,nee) = sum_ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = a*rk_sum(n,nec) + dtsp(n)*residual(n,nec)
      sum_mx = a*rk_sum(n,nmx) + dtsp(n)*residual(n,nmx)
      sum_my = a*rk_sum(n,nmy) + dtsp(n)*residual(n,nmy)
      sum_ee = a*rk_sum(n,nee) + dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nmy) = usp(n,nmy) + b*sum_my
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nmy) = sum_my
      rk_sum(n,nee) = sum_ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,rk_sum,dtsp,residual)
    do m = 1,nq
      do n = 1,np
        !
        sum_eq = a*rk_sum(n,m) + dtsp(n)*residual(n,m)
        !
        usp(n,m) = usp(n,m) + b*sum_eq
        !
        rk_sum(n,m) = sum_eq
        !
      end do
    end do
  end if
  !
end subroutine update_CK_RK4_S4
!
!###############################################################################
!
pure subroutine update_CK_RK4_S5(np,dtsp,residual,usp,rk_sum)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rk_sum(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  real(wp) :: sum_ec,sum_mx,sum_my,sum_mz,sum_ee,sum_eq
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_CK_RK4_S5"
  !
  real(wp), parameter :: a = -1275806237668.0_wp / 842570457699.0_wp
  real(wp), parameter :: b = 2277821191437.0_wp / 14882151754819.0_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = a*rk_sum(n,nec) + dtsp(n)*residual(n,nec)
      sum_mx = a*rk_sum(n,nmx) + dtsp(n)*residual(n,nmx)
      sum_my = a*rk_sum(n,nmy) + dtsp(n)*residual(n,nmy)
      sum_mz = a*rk_sum(n,nmz) + dtsp(n)*residual(n,nmz)
      sum_ee = a*rk_sum(n,nee) + dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nmy) = usp(n,nmy) + b*sum_my
      usp(n,nmz) = usp(n,nmz) + b*sum_mz
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nmy) = sum_my
      rk_sum(n,nmz) = sum_mz
      rk_sum(n,nee) = sum_ee
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,rk_sum,dtsp,residual)
    do n = 1,np
      !
      sum_ec = a*rk_sum(n,nec) + dtsp(n)*residual(n,nec)
      sum_mx = a*rk_sum(n,nmx) + dtsp(n)*residual(n,nmx)
      sum_my = a*rk_sum(n,nmy) + dtsp(n)*residual(n,nmy)
      sum_ee = a*rk_sum(n,nee) + dtsp(n)*residual(n,nee)
      !
      usp(n,nec) = usp(n,nec) + b*sum_ec
      usp(n,nmx) = usp(n,nmx) + b*sum_mx
      usp(n,nmy) = usp(n,nmy) + b*sum_my
      usp(n,nee) = usp(n,nee) + b*sum_ee
      !
      rk_sum(n,nec) = sum_ec
      rk_sum(n,nmx) = sum_mx
      rk_sum(n,nmy) = sum_my
      rk_sum(n,nee) = sum_ee
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,rk_sum,dtsp,residual)
    do m = 1,nq
      do n = 1,np
        !
        sum_eq = a*rk_sum(n,m) + dtsp(n)*residual(n,m)
        !
        usp(n,m) = usp(n,m) + b*sum_eq
        !
        rk_sum(n,m) = sum_eq
        !
      end do
    end do
  end if
  !
end subroutine update_CK_RK4_S5
!
!###############################################################################
!
subroutine update_TVD_RK2(nst,Apply_Filter,Apply_Limiter)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use flowvar, only : residual,usp,uoldsp,dtsp
  use rhs_mod, only : stabilize_usp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK2"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nst == 1) then
    call update_TVD_RK2_S1(n_solpts,dtsp,residual,usp)
  else
    call update_TVD_RK2_S2(n_solpts,dtsp,uoldsp,residual,usp)
  end if
  !
  ! Apply filtering or limiting to stabilize the solution, and then update
  ! the thermodynamic and transport properties at the solution points
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Limiter,caller_name=pname)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine update_TVD_RK2
!
!###############################################################################
!
pure subroutine update_TVD_RK2_S1(np,dtsp,residual,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK2_S1"
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,usp)
    do n = 1,np
      usp(n,nec) = usp(n,nec) + dtsp(n)*residual(n,nec)
      usp(n,nmx) = usp(n,nmx) + dtsp(n)*residual(n,nmx)
      usp(n,nmy) = usp(n,nmy) + dtsp(n)*residual(n,nmy)
      usp(n,nmz) = usp(n,nmz) + dtsp(n)*residual(n,nmz)
      usp(n,nee) = usp(n,nee) + dtsp(n)*residual(n,nee)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,usp)
    do n = 1,np
      usp(n,nec) = usp(n,nec) + dtsp(n)*residual(n,nec)
      usp(n,nmx) = usp(n,nmx) + dtsp(n)*residual(n,nmx)
      usp(n,nmy) = usp(n,nmy) + dtsp(n)*residual(n,nmy)
      usp(n,nee) = usp(n,nee) + dtsp(n)*residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(residual,dtsp,usp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = usp(n,m) + dtsp(n)*residual(n,m)
      end do
    end do
  end if
  !
end subroutine update_TVD_RK2_S1
!
!###############################################################################
!
pure subroutine update_TVD_RK2_S2(np,dtsp,uoldsp,residual,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK2_S2"
  !
  real(wp), parameter :: a = half
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,np
      usp(n,nec) = a * ( uoldsp(n,nec) + usp(n,nec) + dtsp(n)*residual(n,nec) )
      usp(n,nmx) = a * ( uoldsp(n,nmx) + usp(n,nmx) + dtsp(n)*residual(n,nmx) )
      usp(n,nmy) = a * ( uoldsp(n,nmy) + usp(n,nmy) + dtsp(n)*residual(n,nmy) )
      usp(n,nmz) = a * ( uoldsp(n,nmz) + usp(n,nmz) + dtsp(n)*residual(n,nmz) )
      usp(n,nee) = a * ( uoldsp(n,nee) + usp(n,nee) + dtsp(n)*residual(n,nee) )
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,np
      usp(n,nec) = a * ( uoldsp(n,nec) + usp(n,nec) + dtsp(n)*residual(n,nec) )
      usp(n,nmx) = a * ( uoldsp(n,nmx) + usp(n,nmx) + dtsp(n)*residual(n,nmx) )
      usp(n,nmy) = a * ( uoldsp(n,nmy) + usp(n,nmy) + dtsp(n)*residual(n,nmy) )
      usp(n,nee) = a * ( uoldsp(n,nee) + usp(n,nee) + dtsp(n)*residual(n,nee) )
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(residual,dtsp,uoldsp,usp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = a * ( uoldsp(n,m) + usp(n,m) + dtsp(n)*residual(n,m) )
      end do
    end do
  end if
  !
end subroutine update_TVD_RK2_S2
!
!###############################################################################
!
subroutine update_TVD_RK3(nst,Apply_Filter,Apply_Limiter,iter)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use flowvar, only : residual,usp,uoldsp,dtsp
  use rhs_mod, only : stabilize_usp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst,iter
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Scalars ..
  integer :: n
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK3"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nst == 1) then
    call update_TVD_RK3_S1(n_solpts,dtsp,residual,usp)
  else if (nst == 2) then
    call update_TVD_RK3_S2(n_solpts,dtsp,uoldsp,residual,usp)
  else if (nst == 3) then
    call update_TVD_RK3_S3(n_solpts,dtsp,uoldsp,residual,usp)
  end if
  !
  ! Apply filtering or limiting to stabilize the solution, and then update
  ! the thermodynamic and transport properties at the solution points
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Limiter,caller_name=pname)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine update_TVD_RK3
!
!###############################################################################
!
pure subroutine update_TVD_RK3_S1(np,dtsp,residual,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK3_S1"
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,usp)
    do n = 1,np
      usp(n,nec) = usp(n,nec) + dtsp(n)*residual(n,nec)
      usp(n,nmx) = usp(n,nmx) + dtsp(n)*residual(n,nmx)
      usp(n,nmy) = usp(n,nmy) + dtsp(n)*residual(n,nmy)
      usp(n,nmz) = usp(n,nmz) + dtsp(n)*residual(n,nmz)
      usp(n,nee) = usp(n,nee) + dtsp(n)*residual(n,nee)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,usp)
    do n = 1,np
      usp(n,nec) = usp(n,nec) + dtsp(n)*residual(n,nec)
      usp(n,nmx) = usp(n,nmx) + dtsp(n)*residual(n,nmx)
      usp(n,nmy) = usp(n,nmy) + dtsp(n)*residual(n,nmy)
      usp(n,nee) = usp(n,nee) + dtsp(n)*residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(residual,dtsp,usp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = usp(n,m) + dtsp(n)*residual(n,m)
      end do
    end do
  end if
  !
end subroutine update_TVD_RK3_S1
!
!###############################################################################
!
pure subroutine update_TVD_RK3_S2(np,dtsp,uoldsp,residual,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK3_S2"
  !
  real(wp), parameter :: a = three4
  real(wp), parameter :: b = one4
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,np
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + b*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + b*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + b*dtsp(n)*residual(n,nmy)
      usp(n,nmz) = a*uoldsp(n,nmz) + b*usp(n,nmz) + b*dtsp(n)*residual(n,nmz)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + b*dtsp(n)*residual(n,nee)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,np
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + b*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + b*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + b*dtsp(n)*residual(n,nmy)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + b*dtsp(n)*residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(residual,dtsp,uoldsp,usp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = a*uoldsp(n,m) + b*usp(n,m) + b*dtsp(n)*residual(n,m)
      end do
    end do
  end if
  !
end subroutine update_TVD_RK3_S2
!
!###############################################################################
!
pure subroutine update_TVD_RK3_S3(np,dtsp,uoldsp,residual,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK3_S3"
  !
  real(wp), parameter :: a = one3
  real(wp), parameter :: b = two3
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,np
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + b*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + b*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + b*dtsp(n)*residual(n,nmy)
      usp(n,nmz) = a*uoldsp(n,nmz) + b*usp(n,nmz) + b*dtsp(n)*residual(n,nmz)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + b*dtsp(n)*residual(n,nee)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,dtsp,uoldsp,usp)
    do n = 1,np
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + b*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + b*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + b*dtsp(n)*residual(n,nmy)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + b*dtsp(n)*residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(residual,dtsp,uoldsp,usp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = a*uoldsp(n,m) + b*usp(n,m) + b*dtsp(n)*residual(n,m)
      end do
    end do
  end if
  !
end subroutine update_TVD_RK3_S3
!
!###############################################################################
!
subroutine update_TVD_RK4(nst,Apply_Filter,Apply_Limiter)
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use flowvar, only : residual,usp,uoldsp,rk_sum,dtsp
  use rhs_mod, only : stabilize_usp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK4"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  if (nst == 1) then
    call update_TVD_RK4_S1(n_solpts,dtsp,uoldsp,residual,usp)
  else if (nst == 2) then
    call update_TVD_RK4_S2(n_solpts,dtsp,uoldsp,residual,usp)
  else if (nst == 3) then
    call update_TVD_RK4_S3(n_solpts,dtsp,uoldsp,residual,usp,rk_sum)
  else if (nst == 4) then
    call update_TVD_RK4_S4(n_solpts,dtsp,uoldsp,residual,usp,rk_sum)
  else if (nst == 5) then
    call update_TVD_RK4_S5(n_solpts,dtsp,uoldsp,residual,usp,rk_sum)
  end if
  !
  ! Apply filtering or limiting to stabilize the solution, and then update
  ! the thermodynamic and transport properties at the solution points
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Limiter,caller_name=pname)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine update_TVD_RK4
!
!###############################################################################
!
pure subroutine update_TVD_RK4_S1(np,dtsp,uoldsp,residual,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK4_S1"
  !
  real(wp), parameter :: c = 0.39175222700395_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = uoldsp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = uoldsp(n,nmy) + c*dtsp(n)*residual(n,nmy)
      usp(n,nmz) = uoldsp(n,nmz) + c*dtsp(n)*residual(n,nmz)
      usp(n,nee) = uoldsp(n,nee) + c*dtsp(n)*residual(n,nee)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = uoldsp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = uoldsp(n,nmy) + c*dtsp(n)*residual(n,nmy)
      usp(n,nee) = uoldsp(n,nee) + c*dtsp(n)*residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,dtsp,residual)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = uoldsp(n,m) + c*dtsp(n)*residual(n,m)
      end do
    end do
  end if
  !
end subroutine update_TVD_RK4_S1
!
!###############################################################################
!
pure subroutine update_TVD_RK4_S2(np,dtsp,uoldsp,residual,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK4_S2"
  !
  real(wp), parameter :: a = 0.44437049406734_wp
  real(wp), parameter :: b = 0.55562950593266_wp
  real(wp), parameter :: c  = 0.36841059262959_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual)
    do n = 1,np
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + c*dtsp(n)*residual(n,nmy)
      usp(n,nmz) = a*uoldsp(n,nmz) + b*usp(n,nmz) + c*dtsp(n)*residual(n,nmz)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*dtsp(n)*residual(n,nee)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual)
    do n = 1,np
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + c*dtsp(n)*residual(n,nmy)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*dtsp(n)*residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,dtsp,residual)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = a*uoldsp(n,m) + b*usp(n,m) + c*dtsp(n)*residual(n,m)
      end do
    end do
  end if
  !
end subroutine update_TVD_RK4_S2
!
!###############################################################################
!
pure subroutine update_TVD_RK4_S3(np,dtsp,uoldsp,residual,usp,rk_sum)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rk_sum(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  real(wp) :: stage2_ec,stage2_mx,stage2_my,stage2_mz,stage2_ee,stage2_eq
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK4_S3"
  !
  real(wp), parameter :: a = 0.62010185138540_wp
  real(wp), parameter :: b = 0.37989814861460_wp
  real(wp), parameter :: c = 0.25189177424738_wp
  real(wp), parameter :: d = 0.51723167208978_wp
  real(wp), parameter :: e = 0.12759831133288_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual,rk_sum)
    do n = 1,np
      !
      stage2_ec = usp(n,nec)
      stage2_mx = usp(n,nmx)
      stage2_my = usp(n,nmy)
      stage2_mz = usp(n,nmz)
      stage2_ee = usp(n,nee)
      !
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + c*dtsp(n)*residual(n,nmy)
      usp(n,nmz) = a*uoldsp(n,nmz) + b*usp(n,nmz) + c*dtsp(n)*residual(n,nmz)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*dtsp(n)*residual(n,nee)
      !
      rk_sum(n,nec) = d*stage2_ec + e*usp(n,nec)
      rk_sum(n,nmx) = d*stage2_mx + e*usp(n,nmx)
      rk_sum(n,nmy) = d*stage2_my + e*usp(n,nmy)
      rk_sum(n,nmz) = d*stage2_mz + e*usp(n,nmz)
      rk_sum(n,nee) = d*stage2_ee + e*usp(n,nee)
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual,rk_sum)
    do n = 1,np
      !
      stage2_ec = usp(n,nec)
      stage2_mx = usp(n,nmx)
      stage2_my = usp(n,nmy)
      stage2_ee = usp(n,nee)
      !
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + c*dtsp(n)*residual(n,nmy)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*dtsp(n)*residual(n,nee)
      !
      rk_sum(n,nec) = d*stage2_ec + e*usp(n,nec)
      rk_sum(n,nmx) = d*stage2_mx + e*usp(n,nmx)
      rk_sum(n,nmy) = d*stage2_my + e*usp(n,nmy)
      rk_sum(n,nee) = d*stage2_ee + e*usp(n,nee)
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,dtsp,residual,rk_sum)
    do m = 1,nq
      do n = 1,np
        !
        stage2_eq = usp(n,m)
        !
        usp(n,m) = a*uoldsp(n,m) + b*usp(n,m) + c*dtsp(n)*residual(n,m)
        !
        rk_sum(n,m) = d*stage2_eq + e*usp(n,m)
        !
      end do
    end do
  end if
  !
end subroutine update_TVD_RK4_S3
!
!###############################################################################
!
pure subroutine update_TVD_RK4_S4(np,dtsp,uoldsp,residual,usp,rk_sum)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: residual(:,:)
  real(wp), intent(inout) :: usp(:,:)
  real(wp), intent(inout) :: rk_sum(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK4_S4"
  !
  real(wp), parameter :: a = 0.17807995410773_wp
  real(wp), parameter :: b = 0.82192004589227_wp
  real(wp), parameter :: c = 0.54497475021237_wp
  real(wp), parameter :: d = 0.08460416338212_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual,rk_sum)
    do n = 1,np
      !
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + c*dtsp(n)*residual(n,nmy)
      usp(n,nmz) = a*uoldsp(n,nmz) + b*usp(n,nmz) + c*dtsp(n)*residual(n,nmz)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*dtsp(n)*residual(n,nee)
      !
      rk_sum(n,nec) = rk_sum(n,nec) + d*dtsp(n)*residual(n,nec)
      rk_sum(n,nmx) = rk_sum(n,nmx) + d*dtsp(n)*residual(n,nmx)
      rk_sum(n,nmy) = rk_sum(n,nmy) + d*dtsp(n)*residual(n,nmy)
      rk_sum(n,nmz) = rk_sum(n,nmz) + d*dtsp(n)*residual(n,nmz)
      rk_sum(n,nee) = rk_sum(n,nee) + d*dtsp(n)*residual(n,nee)
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,residual,rk_sum)
    do n = 1,np
      !
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*dtsp(n)*residual(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*dtsp(n)*residual(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + c*dtsp(n)*residual(n,nmy)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*dtsp(n)*residual(n,nee)
      !
      rk_sum(n,nec) = rk_sum(n,nec) + d*dtsp(n)*residual(n,nec)
      rk_sum(n,nmx) = rk_sum(n,nmx) + d*dtsp(n)*residual(n,nmx)
      rk_sum(n,nmy) = rk_sum(n,nmy) + d*dtsp(n)*residual(n,nmy)
      rk_sum(n,nee) = rk_sum(n,nee) + d*dtsp(n)*residual(n,nee)
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,dtsp,residual,rk_sum)
    do m = 1,nq
      do n = 1,np
        !
        usp(n,m) = a*uoldsp(n,m) + b*usp(n,m) + c*dtsp(n)*residual(n,m)
        !
        rk_sum(n,m) = rk_sum(n,m) + d*dtsp(n)*residual(n,m)
        !
      end do
    end do
  end if
  !
end subroutine update_TVD_RK4_S4
!
!###############################################################################
!
pure subroutine update_TVD_RK4_S5(np,dtsp,uoldsp,res,usp,rk_sum)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: dtsp(:)
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: res(:,:)
  real(wp),    intent(in) :: rk_sum(:,:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  real(wp) :: c
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_TVD_RK4_S5"
  !
  real(wp), parameter :: a = 0.00683325884039_wp
  real(wp), parameter :: b = 0.34833675773694_wp
  real(wp), parameter :: d = 0.22600748319395_wp
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,res,rk_sum)
    do n = 1,np
      !
      c = d*dtsp(n)
      !
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*res(n,nec) + rk_sum(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*res(n,nmx) + rk_sum(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + c*res(n,nmy) + rk_sum(n,nmy)
      usp(n,nmz) = a*uoldsp(n,nmz) + b*usp(n,nmz) + c*res(n,nmz) + rk_sum(n,nmz)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*res(n,nee) + rk_sum(n,nee)
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,dtsp,res,rk_sum)
    do n = 1,np
      !
      c = d*dtsp(n)
      !
      usp(n,nec) = a*uoldsp(n,nec) + b*usp(n,nec) + c*res(n,nec) + rk_sum(n,nec)
      usp(n,nmx) = a*uoldsp(n,nmx) + b*usp(n,nmx) + c*res(n,nmx) + rk_sum(n,nmx)
      usp(n,nmy) = a*uoldsp(n,nmy) + b*usp(n,nmy) + c*res(n,nmy) + rk_sum(n,nmy)
      usp(n,nee) = a*uoldsp(n,nee) + b*usp(n,nee) + c*res(n,nee) + rk_sum(n,nee)
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,dtsp,res,rk_sum)
    do m = 1,nq
      do n = 1,np
        !
        c = d*dtsp(n)
        !
        usp(n,m) = a*uoldsp(n,m) + b*usp(n,m) + c*res(n,m) + rk_sum(n,m)
        !
      end do
    end do
  end if
  !
end subroutine update_TVD_RK4_S5
!
!###############################################################################
!
subroutine update_Classic_RK(nst,Apply_Filter,Apply_Limiter)
  !
  ! First, store d(utd)/dt for runge kutta time
  ! stepping then update usp
  !
  !.. Use Statements ..
  use geovar,  only : n_solpts
  use ovar,    only : num_rk_stages
  use flowvar, only : rssprk,residual,usp,uoldsp,dtsp
  use rhs_mod, only : stabilize_usp
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: nst
  logical(lk), intent(in) :: Apply_Filter
  logical(lk), intent(in) :: Apply_Limiter
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
  !.. Local Parameters ..
  character(len=*), parameter :: pname = "update_Classic_RK"
  !
continue
  !
  call debug_timer(entering_procedure,pname)
  !
  nr = nme-nmb+1
  !
  ! Store the residual for this Runge-Kutta step in rssprk
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,rssprk)
    do n = 1,n_solpts
      rssprk(n,nec,nst) = residual(n,nec)
      rssprk(n,nmx,nst) = residual(n,nmx)
      rssprk(n,nmy,nst) = residual(n,nmy)
      rssprk(n,nmz,nst) = residual(n,nmz)
      rssprk(n,nee,nst) = residual(n,nee)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(residual,rssprk)
    do n = 1,n_solpts
      rssprk(n,nec,nst) = residual(n,nec)
      rssprk(n,nmx,nst) = residual(n,nmx)
      rssprk(n,nmy,nst) = residual(n,nmy)
      rssprk(n,nee,nst) = residual(n,nee)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(residual,rssprk)
    do m = 1,nq
      do n = 1,n_solpts
        rssprk(n,m,nst) = residual(n,m)
      end do
    end do
  end if
  !
  ! Update usp
  !
  if (num_rk_stages == 1) then
    call update_Classic_RK1_S1(n_solpts,uoldsp,rssprk,dtsp,usp)
  else if (num_rk_stages == 2) then
    if (nst == 1) then
      call update_Classic_RK2_S1(n_solpts,uoldsp,rssprk,dtsp,usp)
    else if (nst == 2) then
      call update_Classic_RK2_S2(n_solpts,uoldsp,rssprk,dtsp,usp)
    end if
  else if (num_rk_stages == 3) then
    if (nst == 1) then
      call update_Classic_RK3_S1(n_solpts,uoldsp,rssprk,dtsp,usp)
    else if (nst == 2) then
      call update_Classic_RK3_S2(n_solpts,uoldsp,rssprk,dtsp,usp)
    else if (nst == 3) then
      call update_Classic_RK3_S3(n_solpts,uoldsp,rssprk,dtsp,usp)
    end if
  else if (num_rk_stages == 3) then
    if (nst == 1) then
      call update_Classic_RK4_S1(n_solpts,uoldsp,rssprk,dtsp,usp)
    else if (nst == 2) then
      call update_Classic_RK4_S2(n_solpts,uoldsp,rssprk,dtsp,usp)
    else if (nst == 3) then
      call update_Classic_RK4_S3(n_solpts,uoldsp,rssprk,dtsp,usp)
    else if (nst == 4) then
      call update_Classic_RK4_S4(n_solpts,uoldsp,rssprk,dtsp,usp)
    end if
  end if
  !
  ! Apply filtering or limiting to stabilize the solution, and then update
  ! the thermodynamic and transport properties at the solution points
  !
  call stabilize_usp(nst,Apply_Filter,Apply_Limiter,caller_name=pname)
  !
  call debug_timer(leaving_procedure,pname)
  !
end subroutine update_Classic_RK
!
!###############################################################################
!
pure subroutine update_Classic_RK1_S1(np,uoldsp,rssprk,dtsp,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp),    intent(in) :: dtsp(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + dtsp(n)*rssprk(n,nmy,1)
      usp(n,nmz) = uoldsp(n,nmz) + dtsp(n)*rssprk(n,nmz,1)
      usp(n,nee) = uoldsp(n,nee) + dtsp(n)*rssprk(n,nee,1)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + dtsp(n)*rssprk(n,nmy,1)
      usp(n,nee) = uoldsp(n,nee) + dtsp(n)*rssprk(n,nee,1)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = uoldsp(n,m) + dtsp(n)*rssprk(n,m,1)
      end do
    end do
  end if
  !
end subroutine update_Classic_RK1_S1
!
!###############################################################################
!
pure subroutine update_Classic_RK2_S1(np,uoldsp,rssprk,dtsp,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp),    intent(in) :: dtsp(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + half*dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + half*dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + half*dtsp(n)*rssprk(n,nmy,1)
      usp(n,nmz) = uoldsp(n,nmz) + half*dtsp(n)*rssprk(n,nmz,1)
      usp(n,nee) = uoldsp(n,nee) + half*dtsp(n)*rssprk(n,nee,1)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + half*dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + half*dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + half*dtsp(n)*rssprk(n,nmy,1)
      usp(n,nee) = uoldsp(n,nee) + half*dtsp(n)*rssprk(n,nee,1)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = uoldsp(n,m) + half*dtsp(n)*rssprk(n,m,1)
      end do
    end do
  end if
  !
end subroutine update_Classic_RK2_S1
!
!###############################################################################
!
pure subroutine update_Classic_RK2_S2(np,uoldsp,rssprk,dtsp,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp),    intent(in) :: dtsp(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + dtsp(n)*rssprk(n,nec,2)
      usp(n,nmx) = uoldsp(n,nmx) + dtsp(n)*rssprk(n,nmx,2)
      usp(n,nmy) = uoldsp(n,nmy) + dtsp(n)*rssprk(n,nmy,2)
      usp(n,nmz) = uoldsp(n,nmz) + dtsp(n)*rssprk(n,nmz,2)
      usp(n,nee) = uoldsp(n,nee) + dtsp(n)*rssprk(n,nee,2)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + dtsp(n)*rssprk(n,nec,2)
      usp(n,nmx) = uoldsp(n,nmx) + dtsp(n)*rssprk(n,nmx,2)
      usp(n,nmy) = uoldsp(n,nmy) + dtsp(n)*rssprk(n,nmy,2)
      usp(n,nee) = uoldsp(n,nee) + dtsp(n)*rssprk(n,nee,2)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = uoldsp(n,m) + dtsp(n)*rssprk(n,m,2)
      end do
    end do
  end if
  !
end subroutine update_Classic_RK2_S2
!
!###############################################################################
!
pure subroutine update_Classic_RK3_S1(np,uoldsp,rssprk,dtsp,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp),    intent(in) :: dtsp(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + one3*dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + one3*dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + one3*dtsp(n)*rssprk(n,nmy,1)
      usp(n,nmz) = uoldsp(n,nmz) + one3*dtsp(n)*rssprk(n,nmz,1)
      usp(n,nee) = uoldsp(n,nee) + one3*dtsp(n)*rssprk(n,nee,1)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + one3*dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + one3*dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + one3*dtsp(n)*rssprk(n,nmy,1)
      usp(n,nee) = uoldsp(n,nee) + one3*dtsp(n)*rssprk(n,nee,1)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = uoldsp(n,m) + one3*dtsp(n)*rssprk(n,m,1)
      end do
    end do
  end if
  !
end subroutine update_Classic_RK3_S1
!
!###############################################################################
!
pure subroutine update_Classic_RK3_S2(np,uoldsp,rssprk,dtsp,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp),    intent(in) :: dtsp(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + two3*dtsp(n)*rssprk(n,nec,2)
      usp(n,nmx) = uoldsp(n,nmx) + two3*dtsp(n)*rssprk(n,nmx,2)
      usp(n,nmy) = uoldsp(n,nmy) + two3*dtsp(n)*rssprk(n,nmy,2)
      usp(n,nmz) = uoldsp(n,nmz) + two3*dtsp(n)*rssprk(n,nmz,2)
      usp(n,nee) = uoldsp(n,nee) + two3*dtsp(n)*rssprk(n,nee,2)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + two3*dtsp(n)*rssprk(n,nec,2)
      usp(n,nmx) = uoldsp(n,nmx) + two3*dtsp(n)*rssprk(n,nmx,2)
      usp(n,nmy) = uoldsp(n,nmy) + two3*dtsp(n)*rssprk(n,nmy,2)
      usp(n,nee) = uoldsp(n,nee) + two3*dtsp(n)*rssprk(n,nee,2)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = uoldsp(n,m) + two3*dtsp(n)*rssprk(n,m,2)
      end do
    end do
  end if
  !
end subroutine update_Classic_RK3_S2
!
!###############################################################################
!
pure subroutine update_Classic_RK3_S3(np,uoldsp,rssprk,dtsp,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp),    intent(in) :: dtsp(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  real(wp) :: c1,c2
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      !
      c1 = one4*dtsp(n)
      c2 = three4*dtsp(n)
      !
      usp(n,nec) = uoldsp(n,nec) + c1*rssprk(n,nec,1) + c2*rssprk(n,nec,3)
      usp(n,nmx) = uoldsp(n,nmx) + c1*rssprk(n,nmx,1) + c2*rssprk(n,nmx,3)
      usp(n,nmy) = uoldsp(n,nmy) + c1*rssprk(n,nmy,1) + c2*rssprk(n,nmy,3)
      usp(n,nmz) = uoldsp(n,nmz) + c1*rssprk(n,nmz,1) + c2*rssprk(n,nmz,3)
      usp(n,nee) = uoldsp(n,nee) + c1*rssprk(n,nee,1) + c2*rssprk(n,nee,3)
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      !
      c1 = one4*dtsp(n)
      c2 = three4*dtsp(n)
      !
      usp(n,nec) = uoldsp(n,nec) + c1*rssprk(n,nec,1) + c2*rssprk(n,nec,3)
      usp(n,nmx) = uoldsp(n,nmx) + c1*rssprk(n,nmx,1) + c2*rssprk(n,nmx,3)
      usp(n,nmy) = uoldsp(n,nmy) + c1*rssprk(n,nmy,1) + c2*rssprk(n,nmy,3)
      usp(n,nee) = uoldsp(n,nee) + c1*rssprk(n,nee,1) + c2*rssprk(n,nee,3)
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        c1 = one4*dtsp(n)
        c2 = three4*dtsp(n)
        !
        usp(n,m) = uoldsp(n,m) + c1*rssprk(n,m,1) + c2*rssprk(n,m,3)
        !
      end do
    end do
  end if
  !
end subroutine update_Classic_RK3_S3
!
!###############################################################################
!
pure subroutine update_Classic_RK4_S1(np,uoldsp,rssprk,dtsp,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp),    intent(in) :: dtsp(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + half*dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + half*dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + half*dtsp(n)*rssprk(n,nmy,1)
      usp(n,nmz) = uoldsp(n,nmz) + half*dtsp(n)*rssprk(n,nmz,1)
      usp(n,nee) = uoldsp(n,nee) + half*dtsp(n)*rssprk(n,nee,1)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + half*dtsp(n)*rssprk(n,nec,1)
      usp(n,nmx) = uoldsp(n,nmx) + half*dtsp(n)*rssprk(n,nmx,1)
      usp(n,nmy) = uoldsp(n,nmy) + half*dtsp(n)*rssprk(n,nmy,1)
      usp(n,nee) = uoldsp(n,nee) + half*dtsp(n)*rssprk(n,nee,1)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = uoldsp(n,m) + half*dtsp(n)*rssprk(n,m,1)
      end do
    end do
  end if
  !
end subroutine update_Classic_RK4_S1
!
!###############################################################################
!
pure subroutine update_Classic_RK4_S2(np,uoldsp,rssprk,dtsp,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp),    intent(in) :: dtsp(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + half*dtsp(n)*rssprk(n,nec,2)
      usp(n,nmx) = uoldsp(n,nmx) + half*dtsp(n)*rssprk(n,nmx,2)
      usp(n,nmy) = uoldsp(n,nmy) + half*dtsp(n)*rssprk(n,nmy,2)
      usp(n,nmz) = uoldsp(n,nmz) + half*dtsp(n)*rssprk(n,nmz,2)
      usp(n,nee) = uoldsp(n,nee) + half*dtsp(n)*rssprk(n,nee,2)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + half*dtsp(n)*rssprk(n,nec,2)
      usp(n,nmx) = uoldsp(n,nmx) + half*dtsp(n)*rssprk(n,nmx,2)
      usp(n,nmy) = uoldsp(n,nmy) + half*dtsp(n)*rssprk(n,nmy,2)
      usp(n,nee) = uoldsp(n,nee) + half*dtsp(n)*rssprk(n,nee,2)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = uoldsp(n,m) + half*dtsp(n)*rssprk(n,m,2)
      end do
    end do
  end if
  !
end subroutine update_Classic_RK4_S2
!
!###############################################################################
!
pure subroutine update_Classic_RK4_S3(np,uoldsp,rssprk,dtsp,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp),    intent(in) :: dtsp(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer :: m,n,nr
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + dtsp(n)*rssprk(n,nec,3)
      usp(n,nmx) = uoldsp(n,nmx) + dtsp(n)*rssprk(n,nmx,3)
      usp(n,nmy) = uoldsp(n,nmy) + dtsp(n)*rssprk(n,nmy,3)
      usp(n,nmz) = uoldsp(n,nmz) + dtsp(n)*rssprk(n,nmz,3)
      usp(n,nee) = uoldsp(n,nee) + dtsp(n)*rssprk(n,nee,3)
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      usp(n,nec) = uoldsp(n,nec) + dtsp(n)*rssprk(n,nec,3)
      usp(n,nmx) = uoldsp(n,nmx) + dtsp(n)*rssprk(n,nmx,3)
      usp(n,nmy) = uoldsp(n,nmy) + dtsp(n)*rssprk(n,nmy,3)
      usp(n,nee) = uoldsp(n,nee) + dtsp(n)*rssprk(n,nee,3)
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp)
    do m = 1,nq
      do n = 1,np
        usp(n,m) = uoldsp(n,m) + dtsp(n)*rssprk(n,m,3)
      end do
    end do
  end if
  !
end subroutine update_Classic_RK4_S3
!
!###############################################################################
!
pure subroutine update_Classic_RK4_S4(np,uoldsp,rssprk,dtsp,usp)
  !
  !.. Formal Arguments ..
  integer,     intent(in) :: np
  real(wp),    intent(in) :: uoldsp(:,:)
  real(wp),    intent(in) :: rssprk(:,:,:)
  real(wp),    intent(in) :: dtsp(:)
  real(wp), intent(inout) :: usp(:,:)
  !
  !.. Local Scalars ..
  integer  :: m,n,nr
  real(wp) :: c1,c2
  !
continue
  !
  nr = nme-nmb+1
  !
  if (nr == 3 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      !
      c1 = one6*dtsp(n)
      c2 = one3*dtsp(n)
      !
      usp(n,nec) = uoldsp(n,nec) + c1*rssprk(n,nec,1) + c2*rssprk(n,nec,2) &
                                 + c2*rssprk(n,nec,3) + c1*rssprk(n,nec,4)
      usp(n,nmx) = uoldsp(n,nmx) + c1*rssprk(n,nmx,1) + c2*rssprk(n,nmx,2) &
                                 + c2*rssprk(n,nmx,3) + c1*rssprk(n,nmx,4)
      usp(n,nmy) = uoldsp(n,nmy) + c1*rssprk(n,nmy,1) + c2*rssprk(n,nmy,2) &
                                 + c2*rssprk(n,nmy,3) + c1*rssprk(n,nmy,4)
      usp(n,nmz) = uoldsp(n,nmz) + c1*rssprk(n,nmz,1) + c2*rssprk(n,nmz,2) &
                                 + c2*rssprk(n,nmz,3) + c1*rssprk(n,nmz,4)
      usp(n,nee) = uoldsp(n,nee) + c1*rssprk(n,nee,1) + c2*rssprk(n,nee,2) &
                                 + c2*rssprk(n,nee,3) + c1*rssprk(n,nee,4)
      !
    end do
  else if (nr == 2 .and. ncs == 0 .and. .not.force_nd_subroutines) then
    !$ACC PARALLEL LOOP PRESENT(usp,uoldsp,rssprk,dtsp)
    do n = 1,np
      !
      c1 = one6*dtsp(n)
      c2 = one3*dtsp(n)
      !
      usp(n,nec) = uoldsp(n,nec) + c1*rssprk(n,nec,1) + c2*rssprk(n,nec,2) &
                                 + c2*rssprk(n,nec,3) + c1*rssprk(n,nec,4)
      usp(n,nmx) = uoldsp(n,nmx) + c1*rssprk(n,nmx,1) + c2*rssprk(n,nmx,2) &
                                 + c2*rssprk(n,nmx,3) + c1*rssprk(n,nmx,4)
      usp(n,nmy) = uoldsp(n,nmy) + c1*rssprk(n,nmy,1) + c2*rssprk(n,nmy,2) &
                                 + c2*rssprk(n,nmy,3) + c1*rssprk(n,nmy,4)
      usp(n,nee) = uoldsp(n,nee) + c1*rssprk(n,nee,1) + c2*rssprk(n,nee,2) &
                                 + c2*rssprk(n,nee,3) + c1*rssprk(n,nee,4)
      !
    end do
  else
    !$ACC PARALLEL LOOP COLLAPSE(2) PRESENT(usp,uoldsp,rssprk,dtsp)
    do m = 1,nq
      do n = 1,np
        !
        c1 = one6*dtsp(n)
        c2 = one3*dtsp(n)
        !
        usp(n,m) = uoldsp(n,m) + c1*rssprk(n,m,1) + c2*rssprk(n,m,2) &
                               + c2*rssprk(n,m,3) + c1*rssprk(n,m,4)
        !
      end do
    end do
  end if
  !
end subroutine update_Classic_RK4_S4
!
!###############################################################################
!
#include "Functions/debug_timer.f90"
!
!###############################################################################
!
end module erk_mod
