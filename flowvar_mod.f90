module flowvar
  !
  use module_kind_types, only : vbp
  !
  implicit none
  !
  private
  !
 !include "kind_defs.h"
  !
  !      usp(4,ksp2d,-nfbnd:ncell)
  !      vsp(4,ksp2d,-nfbnd:ncell)
  ! uspstart(4,ksp2d,-nfbnd:ncell)
  ! residual(4,ksp2d,-nfbnd:ncell)
  !   uoldsp(4,ksp2d,-nfbnd:ncell)
  !   VSPbe_(6,ksp1d,nfbnd)
  !   WSPbe_(4,ksp1d,nfbnd)
  !
  real(vbp), public, save, allocatable, dimension(:,:) :: usp
  real(vbp), public, save, allocatable, dimension(:,:) :: vsp
  real(vbp), public, save, allocatable, dimension(:,:) :: uspstart
  real(vbp), public, save, allocatable, dimension(:,:) :: residual
  real(vbp), public, save, allocatable, dimension(:,:) :: uoldsp
  real(vbp), public, save, allocatable, dimension(:,:) :: usp_sum
  real(vbp), public, save, allocatable, dimension(:,:) :: global_usp
  real(vbp), public, save, allocatable, dimension(:,:) :: global_vsp
  !
  ! faceflx(1:neqns,1:2,1:n_sp1de,1:nface) ! this is the normal face flux
  ! faceusp(1:neqns,1:2,1:n_sp1de,1:nface)
  real(vbp), public, save, allocatable, dimension(:,:,:,:) :: faceusp
  real(vbp), public, save, allocatable, dimension(:,:,:,:) :: faceflx
  !
 !real(vbp), public, save, allocatable, dimension(:,:,:) :: VSPbe_
 !real(vbp), public, save, allocatable, dimension(:,:,:) :: WSPbe_
  !
  ! rssprk(4,ksp2d,-nfbnd:ncell,num_stages)
  !
  real(vbp), public, save, allocatable, dimension(:,:,:) :: rssprk
  !
end module flowvar
!
!###############################################################################
!
