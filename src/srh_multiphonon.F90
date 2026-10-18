!!--------------------------------------------------------------------------!
!! libNEGF: a general library for Non-Equilibrium Greens functions.         !
!! Copyright (C) 2012 - 2026                                                !
!!                                                                          !
!! Skeleton SRH / multiphonon defect interaction (Stage 4).                 !
!! Stores defect parameters; Sigma methods are no-ops until physics is      !
!! filled in. Interface matches TInteraction so SCBA can register it.       !
!!--------------------------------------------------------------------------!

module srh_multiphonon

  use ln_precision, only : dp
  use interactions, only : TInteraction
  use ln_inelastic, only : TInelastic
  use ln_structure, only : TStruct_info
  use mat_def, only : c_dns, z_dns

  implicit none
  private

  public :: TSrhMultiphonon
  public :: TSrhMultiphonon_create
  public :: TSrhMultiphonon_init

  !> Localized deep-defect SRH interaction (skeleton).
  !! One defect site for now; Cn/Cp capture coefficients; full multiphonon later.
  type, extends(TInelastic) :: TSrhMultiphonon

    !> Defect energy relative to midgap reference (eV), zero-field
    real(dp) :: Ed0 = 0.0_dp
    !> Volume defect density (cm^-3) or effective density in local cell
    real(dp) :: rho_d = 0.0_dp
    !> Electron / hole capture coefficients (cm^3/s)
    real(dp) :: Cn = 0.0_dp
    real(dp) :: Cp = 0.0_dp
    !> Temperature (K)
    real(dp) :: temperature = 300.0_dp
    !> Principal-layer index of the defect (1-based); 0 = unset / not mapped yet
    integer :: pl_index = 0
    !> Local orbital index within that PL (1-based); 0 = unset
    integer :: orb_index = 0
    !> Steady-state defect occupation (placeholder)
    real(dp) :: fd = 0.5_dp
    !> Last net recombination rate diagnostic (placeholder, volume units)
    real(dp) :: r_net = 0.0_dp

  contains

    procedure :: add_sigma_r_sp
    procedure :: add_sigma_r_dp
    procedure :: add_sigma_n_sp
    procedure :: add_sigma_n_dp
    procedure :: get_sigma_n_blk
    procedure :: get_sigma_n_mat
    procedure :: set_Gr
    procedure :: set_Gn
    procedure :: compute_Sigma_r
    procedure :: compute_Sigma_n
    procedure :: destroy_Sigma_r
    procedure :: destroy_Sigma_n
    procedure :: destroy

  end type TSrhMultiphonon

contains

  subroutine TSrhMultiphonon_create(this)
    class(TInteraction), allocatable :: this
    allocate(TSrhMultiphonon::this)
  end subroutine TSrhMultiphonon_create

  !> Initialize skeleton SRH interaction.
  !! @param struct  device structure (PL partitioning)
  !! @param Ed0     defect energy (eV)
  !! @param rho_d   defect density
  !! @param Cn, Cp  capture coefficients (cm^3/s)
  !! @param temperature  Kelvin
  !! @param pl_index     principal layer index (1-based), 0 if unknown
  !! @param orb_index    orbital within PL (1-based), 0 if unknown
  !! @param niter        SCBA iteration cap
  subroutine TSrhMultiphonon_init(this, struct, Ed0, rho_d, Cn, Cp, temperature, &
      & pl_index, orb_index, niter)
    type(TSrhMultiphonon), intent(inout) :: this
    type(TStruct_info), intent(in) :: struct
    real(dp), intent(in) :: Ed0, rho_d, Cn, Cp, temperature
    integer, intent(in) :: pl_index, orb_index, niter

    this%descriptor = "SRH multiphonon defect interaction (skeleton)"
    this%struct = struct
    this%Ed0 = Ed0
    this%rho_d = rho_d
    this%Cn = Cn
    this%Cp = Cp
    this%temperature = temperature
    this%pl_index = pl_index
    this%orb_index = orb_index
    this%scba_niter = niter
    this%wq = 0.0_dp
    this%fd = 0.5_dp
    this%r_net = 0.0_dp

  end subroutine TSrhMultiphonon_init

  ! ---- deferred TInteraction methods: skeleton no-ops --------------------

  subroutine add_sigma_r_sp(this, esh, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(c_dns), dimension(:,:), intent(inout) :: esh
    integer, intent(in), optional :: en_index, k_index, spin
    if (this%scba_iter .eq. 0) return
    ! skeleton: no Sigma_r contribution yet
  end subroutine add_sigma_r_sp

  subroutine add_sigma_r_dp(this, esh, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(z_dns), dimension(:,:), intent(inout) :: esh
    integer, intent(in), optional :: en_index, k_index, spin
    if (this%scba_iter .eq. 0) return
  end subroutine add_sigma_r_dp

  subroutine add_sigma_n_sp(this, sigma, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(c_dns), dimension(:,:), intent(inout) :: sigma
    integer, intent(in), optional :: en_index, k_index, spin
    if (this%scba_iter .eq. 0) return
  end subroutine add_sigma_n_sp

  subroutine add_sigma_n_dp(this, sigma, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(z_dns), dimension(:,:), intent(inout) :: sigma
    integer, intent(in), optional :: en_index, k_index, spin
    if (this%scba_iter .eq. 0) return
  end subroutine add_sigma_n_dp

  subroutine get_sigma_n_blk(this, blk_sigma_n, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(z_dns), dimension(:,:), intent(inout) :: blk_sigma_n
    integer, intent(in), optional :: en_index, k_index, spin
  end subroutine get_sigma_n_blk

  subroutine get_sigma_n_mat(this, sigma_n, ii, jj, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(z_dns), intent(inout) :: sigma_n
    integer, intent(in) :: ii, jj
    integer, intent(in), optional :: en_index, k_index, spin
  end subroutine get_sigma_n_mat

  subroutine set_Gr(this, Gr, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(z_dns), dimension(:,:), intent(in) :: Gr
    integer, intent(in), optional :: en_index, k_index, spin
    if (this%scba_iter > this%scba_niter) return
    ! skeleton: will read local block at pl_index later
  end subroutine set_Gr

  subroutine set_Gn(this, Gn, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(z_dns), dimension(:,:), intent(in) :: Gn
    integer, intent(in), optional :: en_index, k_index, spin
    if (this%scba_iter > this%scba_niter) return
  end subroutine set_Gn

  subroutine compute_Sigma_r(this, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    integer, intent(in), optional :: en_index, k_index, spin
  end subroutine compute_Sigma_r

  subroutine compute_Sigma_n(this, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    integer, intent(in), optional :: en_index, k_index, spin
  end subroutine compute_Sigma_n

  subroutine destroy_Sigma_r(this)
    class(TSrhMultiphonon) :: this
  end subroutine destroy_Sigma_r

  subroutine destroy_Sigma_n(this)
    class(TSrhMultiphonon) :: this
  end subroutine destroy_Sigma_n

  subroutine destroy(this)
    class(TSrhMultiphonon) :: this
    call this%destroy_Sigma_r()
    call this%destroy_Sigma_n()
  end subroutine destroy

end module srh_multiphonon
