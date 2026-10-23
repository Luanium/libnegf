!!--------------------------------------------------------------------------!
!! libNEGF: a general library for Non-Equilibrium Greens functions.         !
!! Copyright (C) 2012 - 2026                                                !
!!                                                                          !
!! SRH defect-mediated recombination via multiphonon relaxation (NEGF).     !
!!                                                                          !
!! References:                                                              !
!!   [A] U. Aeberhard, PRB 99, 125302 (2019)  Eqs. 45-52, 53-59, App. A, B  !
!!   [P] A. Pilotto et al., Solid-State Electron. 198, 108469 (2022) Eq. 4-7!
!!                                                                          !
!! One interaction object = one defect SPECIES (same level above Ev, same   !
!! density rho_d, same capture coefficients, same phonon ladder) placed on  !
!! one or MANY SITES.  A site is n_orb consecutive matrix rows inside one   !
!! principal layer (pl_index; rows orb_index .. orb_index+n_orb-1).  Every  !
!! site has its own band edges Ev, Ec (they follow the electrostatic        !
!! potential and the alloy composition), hence its own defect energy        !
!! Ed = Ev + Ed_rel, its own calibration, its own occupation fd.            !
!! A uniform defect density in a region = one site on every atom of the     !
!! region (the delta(z-z0) of [A] Eq. 63 smeared over the atom cell).       !
!!                                                                          !
!! Per SCBA iteration and for every site (delta-like defect, [A] Eq. 8):    !
!!  1. local band Green's functions G^n, G^p = A - G^n at the site         !
!!     (average over the n_orb orbitals, then over k) from the libNEGF      !
!!     caches (full G: conduction and valence band are both contained).     !
!!  2. defect self-energies (Eqs. 45-48 of [A]), ladder l = 0..lmax:        !
!!       S^n_dc = sum_l Mn(l)              G^n(E + l hw)                    !
!!       S^p_dc = sum_l Mn(l) e^{-b l hw}  G^p(E + l hw)                    !
!!       S^n_dv = sum_l Mp(l) e^{-b l hw}  G^n(E - l hw)                    !
!!       S^p_dv = sum_l Mp(l)              G^p(E - l hw)                    !
!!     averaged over the (broadened) defect line.                           !
!!  3. defect occupation (Eq. 59 of [A]):                                   !
!!       fd = (S^n_dc + S^n_dv) / (S^n_dc + S^p_dc + S^n_dv + S^p_dv)      !
!!  4. net recombination per defect (Eqs. 55,58; closed form, no cancel.)  !
!!       R = (S^n_dc S^p_dv - S^p_dc S^n_dv) / (Gamma hbar),  U = rho_d R   !
!!  5. band self-energies from the DEFECT Green's function (Eqs. 49-52),    !
!!     G^n_d = fd A_d, G^p_d = (1-fd) A_d, same value on every orbital of   !
!!     the site:                                                           !
!!       S^n_cd(E) = x_d sum_l Mn(l) e^{-b l hw} fd     A_d(E - l hw)       !
!!       S^p_cd(E) = x_d sum_l Mn(l)             (1-fd) A_d(E - l hw)       !
!!       S^n_vd(E) = x_d sum_l Mp(l)             fd     A_d(E + l hw)       !
!!       S^p_vd(E) = x_d sum_l Mp(l) e^{-b l hw} (1-fd) A_d(E + l hw)       !
!!     with x_d = rho_d * cell_volume * weight (defects per orbital cell).  !
!!     S^r(E) = Lambda(E) - (i/2)[S^n(E)+S^p(E)], where Lambda is the      !
!!     principal-value (Kramers-Kronig) transform of the broadening.        !
!!                                                                          !
!! Storage: the band self-energies are DIAGONAL and independent of k, so    !
!! only Sn_b(E), Sp_b(E) per site are stored (no dense PL blocks).          !
!!                                                                          !
!! Conventions: libNEGF G_n = -i G^< (Hermitian, real positive diagonal),   !
!! G_p = i G^> = A - G_n, A = -2 Im G_r.  Sigma_n, Sigma_p likewise.        !
!!                                                                          !
!! Coupling strength: M(l) = hbar kappa shape(l) g_s / (2 pi V_cell), with  !
!! shape(l) the finite-temperature ladder weight of [P] Eq. 6 and kappa     !
!! fixed by the thermal-average calibration of [A] App. B (Eq. B1/B2) with  !
!! analytic parabolic bands: c_n = Cn = sigma_n v_th (c_p = Cp likewise).   !
!! Hole weight uses (l+S)^2/S, electron weight (l-S)^2/S ([A] Eq. A5, A10).  !
!! The calibration is done per site (it depends on the local gap).          !
!!                                                                          !
!! MPI: local k contributions are reduced over the k communicator and the  !
!! resulting energy blocks over the energy communicator before evaluating   !
!! the shifted phonon ladder.  A uniform energy grid is required; the       !
!! defect energy is fixed in absolute energy by local band edges supplied   !
!! when the site is added.                                                   !
!!--------------------------------------------------------------------------!

module srh_multiphonon

  use ln_precision, only : dp
  use ln_messages, only : error_msg
  use mpi_globals, only : id0
  use interactions, only : TInteraction
  use ln_inelastic, only : TInelastic
  use ln_structure, only : TStruct_info
  use mat_def, only : c_dns, z_dns
  use ln_cache, only : TMatrixCacheMem, TMatLabel
#:if defined("MPI")
  ! MPIFX collectives use the legacy integer MPI operation handles, not the
  ! mpi_f08 derived MPI_Op handles.
  use mpi, only : MPI_SUM
  use libmpifx_module, only : mpifx_comm, mpifx_allreduceip
#:endif

  implicit none
  private

  public :: TSrhMultiphonon
  public :: TSrhSite
  public :: TSrhMultiphonon_create
  public :: TSrhMultiphonon_init

  real(dp), parameter :: PI = 3.14159265358979323846_dp
  real(dp), parameter :: HBAR_EV_S = 6.582119569e-16_dp
  real(dp), parameter :: HBAR_J_S = 1.054571817e-34_dp
  real(dp), parameter :: M0_KG = 9.1093837015e-31_dp
  real(dp), parameter :: Q_C = 1.602176634e-19_dp
  real(dp), parameter :: KB_EV = 8.617333262145e-5_dp
  !> spin degeneracy of the band states (defect level itself is spinless, [A] note 35)
  real(dp), parameter :: G_SPIN = 2.0_dp
  !> maximum number of phonons in the ladder
  integer, parameter :: LMAX_RUNGS = 200
  !> defect line is represented within +- LINE_SIGMAS * width
  real(dp), parameter :: LINE_SIGMAS = 6.0_dp

  !> One defect site (n_orb consecutive rows of one principal layer)
  type :: TSrhSite
    integer :: pl = 0
    integer :: orb = 0
    integer :: n_orb = 1
    !> multiplicity (a site standing for `weight` equivalent atoms)
    real(dp) :: weight = 1.0_dp
    !> local band edges, absolute energies on the NEGF axis [eV]
    real(dp) :: Ev = 0.0_dp
    real(dp) :: Ec = 0.0_dp
    real(dp) :: Eg = 0.0_dp
    real(dp) :: Ed_abs = 0.0_dp
    !> M_dc(l), M_dv(l) [eV^2]
    real(dp), allocatable :: Mn(:)
    real(dp), allocatable :: Mp(:)
    !> defect line on the energy grid
    integer :: jlo = 0
    integer :: jhi = 0
    real(dp), allocatable :: aline(:)
    !> state
    real(dp) :: fd = 0.5_dp
    real(dp) :: omfd = 0.5_dp
    real(dp) :: R_def = 0.0_dp
    real(dp) :: U_rec = 0.0_dp
    real(dp) :: Sn_dc = 0.0_dp
    real(dp) :: Sp_dc = 0.0_dp
    real(dp) :: Sn_dv = 0.0_dp
    real(dp) :: Sp_dv = 0.0_dp
    !> band self-energies at the site orbitals as function of energy index [eV]
    real(dp), allocatable :: Sn_b(:)
    real(dp), allocatable :: Sp_b(:)
    !> conduction-ladder part only of Sn_b, Sp_b (valence part = total - this);
    !! needed to split the band-side flux into its CB and VB contributions
    real(dp), allocatable :: Sn_bc(:)
    real(dp), allocatable :: Sp_bc(:)
    !> band-side net electron flux per defect, from the stored band self-energies
    !! and the current G: rb_c = CB -> defect, rb_v = VB -> defect [1/s]
    real(dp) :: rb_c = 0.0_dp
    real(dp) :: rb_v = 0.0_dp
    !> |rb_c - Rc|, |rb_v + Rv| relative to the gross flux (consistency residual)
    real(dp) :: res_band = 0.0_dp
    !> dispersive real part of the retarded band self-energy [eV]
    real(dp), allocatable :: Sr_b(:)
  end type TSrhSite

  type, extends(TInelastic) :: TSrhMultiphonon

    ! ---- species parameters (common to all sites) -------------------------
    !> defect level above the local valence band edge [eV]
    real(dp) :: Ed_rel = 0.0_dp
    !> defect density [cm^-3]
    real(dp) :: rho_d = 0.0_dp
    !> capture coefficients sigma*v_th [cm^3/s]
    real(dp) :: Cn = 0.0_dp
    real(dp) :: Cp = 0.0_dp
    real(dp) :: temperature = 300.0_dp
    !> phonon energy as given by the user [eV]
    real(dp) :: hw_in = 0.0_dp
    !> Huang-Rhys factor
    real(dp) :: S_hr = 0.0_dp
    !> effective masses [m0] (only used for the analytic calibration)
    real(dp) :: mass_c = 0.0_dp
    real(dp) :: mass_v = 0.0_dp
    !> volume of the cell represented by one orbital [cm^3]
    real(dp) :: cell_vol = 0.0_dp
    !> requested defect line width [eV]
    real(dp) :: width_in = 1.0e-3_dp

    ! ---- sites -------------------------------------------------------------
    integer :: n_sites = 0
    type(TSrhSite), allocatable :: sites(:)

    ! ---- derived: energy grid ---------------------------------------------
    logical :: grid_ready = .false.
    real(dp) :: Emin = 0.0_dp
    real(dp) :: dE = 0.0_dp
    integer :: nE = 0
    !> number of energy points retained in the local Green-function cache
    integer :: nE_local = 0
    !> first global energy index represented by the local Green-function cache
    integer :: iE_first = 1
    integer :: iEshift = 0
    !> rounded phonon energy actually used everywhere (iEshift*dE) [eV]
    real(dp) :: sh = 0.0_dp
    real(dp) :: kT = 0.0_dp
    real(dp) :: fB = 0.0_dp
    real(dp) :: width = 0.0_dp
    real(dp), allocatable :: kw(:)
#:if defined("MPI")
    !> communicator over energy blocks and communicator over local k blocks
    type(mpifx_comm) :: energy_comm
    type(mpifx_comm) :: k_comm
#:endif

    ! ---- derived: ladder (common) -----------------------------------------
    !> ladder shapes rescaled to max = 1
    real(dp), allocatable :: shape_e(:)
    real(dp), allocatable :: shape_h(:)
    !> Boltzmann factor exp(-l sh / kT)
    real(dp), allocatable :: bolt(:)

  contains

    procedure :: add_site
    procedure :: get_site_results
    procedure :: get_site_spectral
    procedure :: get_site_retarded_spectral
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
    procedure :: set_EnGrid
#:if defined("MPI")
    procedure :: set_mpi_communicators
#:endif

  end type TSrhMultiphonon

contains

  !--------------------------------------------------------------------------
  subroutine TSrhMultiphonon_create(this)
    class(TInteraction), allocatable :: this
    allocate(TSrhMultiphonon::this)
  end subroutine TSrhMultiphonon_create

  !--------------------------------------------------------------------------
  !> Store the species parameters (sites are added with add_site; the energy
  !! grid is set later by set_EnGrid)
  subroutine TSrhMultiphonon_init(this, struct, Ed_rel, rho_d, Cn, Cp, temperature, &
      & niter, phonon_energy, huang_rhys, mass_c, mass_v, cell_volume, line_width)
    type(TSrhMultiphonon), intent(inout) :: this
    type(TStruct_info), intent(in) :: struct
    real(dp), intent(in) :: Ed_rel, rho_d, Cn, Cp, temperature
    integer, intent(in) :: niter
    real(dp), intent(in) :: phonon_energy, huang_rhys, mass_c, mass_v
    real(dp), intent(in) :: cell_volume, line_width

    if (phonon_energy <= 0.0_dp) call error_msg('SRH defect: phonon_frequency must be > 0 (eV)')
    if (huang_rhys <= 0.0_dp) call error_msg('SRH defect: huang_rhys must be > 0')
    if (mass_c <= 0.0_dp .or. mass_v <= 0.0_dp) &
        & call error_msg('SRH defect: effective_mass_c and effective_mass_v must be > 0 (units of m0)')
    if (Ed_rel <= 0.0_dp) call error_msg('SRH defect: defect_energy (above Ev) must be > 0')
    if (cell_volume <= 0.0_dp) call error_msg('SRH defect: cell volume must be > 0')
    if (rho_d <= 0.0_dp .or. Cn <= 0.0_dp .or. Cp <= 0.0_dp) &
        & call error_msg('SRH defect: defect_density, capture_n, capture_p must be > 0')
    if (temperature <= 0.0_dp) call error_msg('SRH defect: temperature must be > 0')

    this%descriptor = "SRH multiphonon defect (NEGF, Aeberhard 2019)"
    this%struct = struct
    this%Ed_rel = Ed_rel
    this%rho_d = rho_d
    this%Cn = Cn
    this%Cp = Cp
    this%temperature = temperature
    this%scba_niter = niter
    this%hw_in = phonon_energy
    this%S_hr = huang_rhys
    this%mass_c = mass_c
    this%mass_v = mass_v
    this%cell_vol = cell_volume
    if (line_width > 0.0_dp) this%width_in = line_width
    this%kT = KB_EV * temperature

    this%n_sites = 0
    this%grid_ready = .false.
    ! wq is used by the library for the energy-window bookkeeping of inelastic models
    this%wq = phonon_energy

    if (.not. associated(this%sigma_r)) then
      allocate(this%sigma_r, source=TMatrixCacheMem(tagname='SRH_Sigma_r'))
    end if
    if (.not. associated(this%sigma_n)) then
      allocate(this%sigma_n, source=TMatrixCacheMem(tagname='SRH_Sigma_n'))
    end if

  end subroutine TSrhMultiphonon_init

  !--------------------------------------------------------------------------
  !> Append one defect site.  Ev_abs, Ec_abs: local band edges on the NEGF
  !! energy axis [eV]; weight: multiplicity of the site (>= 1, usually 1).
  subroutine add_site(this, pl_index, orb_index, n_orb, Ev_abs, Ec_abs, weight)
    class(TSrhMultiphonon), intent(inout) :: this
    integer, intent(in) :: pl_index, orb_index, n_orb
    real(dp), intent(in) :: Ev_abs, Ec_abs, weight

    type(TSrhSite), allocatable :: tmp(:)
    integer :: n

    if (pl_index < 1 .or. orb_index < 1 .or. n_orb < 1) &
        & call error_msg('SRH defect: site needs pl_index >= 1, orb_index >= 1, n_orb >= 1')
    if (Ec_abs - Ev_abs <= 0.0_dp) &
        & call error_msg('SRH defect: band_edge_c must be above band_edge_v at the defect site')
    if (this%Ed_rel >= Ec_abs - Ev_abs) &
        & call error_msg('SRH defect: defect_energy (above Ev) must lie inside the local gap of every site')
    if (weight <= 0.0_dp) call error_msg('SRH defect: site weight must be > 0')

    n = this%n_sites
    allocate(tmp(n + 1))
    if (n > 0) tmp(1:n) = this%sites(1:n)
    tmp(n + 1)%pl = pl_index
    tmp(n + 1)%orb = orb_index
    tmp(n + 1)%n_orb = n_orb
    tmp(n + 1)%weight = weight
    tmp(n + 1)%Ev = Ev_abs
    tmp(n + 1)%Ec = Ec_abs
    tmp(n + 1)%Eg = Ec_abs - Ev_abs
    tmp(n + 1)%Ed_abs = Ev_abs + this%Ed_rel
    call move_alloc(tmp, this%sites)
    this%n_sites = n + 1
    this%grid_ready = .false.

  end subroutine add_site

  !--------------------------------------------------------------------------
  !> res(1:8) = fd, 1-fd, R_def [1/s], U = rho_d*R_def [cm^-3 s^-1],
  !!            Sn_dc, Sp_dc, Sn_dv, Sp_dv [eV];  res = 0 if isite out of range
  subroutine get_site_results(this, isite, res, ok)
    class(TSrhMultiphonon), intent(in) :: this
    integer, intent(in) :: isite
    real(dp), intent(out) :: res(8)
    logical, intent(out) :: ok

    res = 0.0_dp
    ok = .false.
    if (isite < 1 .or. isite > this%n_sites) return
    res(1) = this%sites(isite)%fd
    res(2) = this%sites(isite)%omfd
    res(3) = this%sites(isite)%R_def
    res(4) = this%sites(isite)%U_rec
    res(5) = this%sites(isite)%Sn_dc
    res(6) = this%sites(isite)%Sp_dc
    res(7) = this%sites(isite)%Sn_dv
    res(8) = this%sites(isite)%Sp_dv
    ok = .true.

  end subroutine get_site_results

  !--------------------------------------------------------------------------
  !> Band self-energies Sn_b(E), Sp_b(E) [eV] on the global energy grid for one site.
  !! nE_out / Emin / dE describe the grid (E_i = Emin + (i-1)*dE, i=1..nE_out).
  !! Sn, Sp must be size >= nE_out (caller allocates).  ok=.false. if site missing
  !! or spectral arrays not yet computed.
  subroutine get_site_spectral(this, isite, nE_out, Emin, dE, Sn, Sp, ok)
    class(TSrhMultiphonon), intent(in) :: this
    integer, intent(in) :: isite
    integer, intent(out) :: nE_out
    real(dp), intent(out) :: Emin, dE
    real(dp), intent(out) :: Sn(*), Sp(*)
    logical, intent(out) :: ok

    integer :: n, i

    nE_out = 0
    Emin = 0.0_dp
    dE = 0.0_dp
    ok = .false.
    if (isite < 1 .or. isite > this%n_sites) return
    if (.not. allocated(this%sites(isite)%Sn_b) .or. .not. allocated(this%sites(isite)%Sp_b)) return
    n = size(this%sites(isite)%Sn_b)
    if (n < 1) return
    nE_out = n
    Emin = this%Emin
    dE = this%dE
    do i = 1, n
      Sn(i) = this%sites(isite)%Sn_b(i)
      Sp(i) = this%sites(isite)%Sp_b(i)
    end do
    ok = .true.

  end subroutine get_site_spectral

  !--------------------------------------------------------------------------
  !> Dispersive real component Sr_b(E) [eV] of the causal retarded band
  !! self-energy for one site.  The complete retarded component is
  !! Sigma^r_b = Sr_b - i (Sn_b + Sp_b)/2.
  subroutine get_site_retarded_spectral(this, isite, nE_out, Emin, dE, Sr, ok)
    class(TSrhMultiphonon), intent(in) :: this
    integer, intent(in) :: isite
    integer, intent(out) :: nE_out
    real(dp), intent(out) :: Emin, dE
    real(dp), intent(out) :: Sr(*)
    logical, intent(out) :: ok

    integer :: n, i

    nE_out = 0
    Emin = 0.0_dp
    dE = 0.0_dp
    ok = .false.
    if (isite < 1 .or. isite > this%n_sites) return
    if (.not. allocated(this%sites(isite)%Sr_b)) return
    n = size(this%sites(isite)%Sr_b)
    if (n < 1) return
    nE_out = n
    Emin = this%Emin
    dE = this%dE
    do i = 1, n
      Sr(i) = this%sites(isite)%Sr_b(i)
    end do
    ok = .true.

  end subroutine get_site_retarded_spectral

  !--------------------------------------------------------------------------
  !> Called by libNEGF before the SCBA loop. Builds everything that depends on
  !! the energy grid: rounded phonon shift, ladder shapes, per-site calibration
  !! and defect line, normalised k weights.
  subroutine set_EnGrid(this, Emin, deltaE, nE_global, nE_local, kweights, local_kindex)
    class(TSrhMultiphonon), intent(inout) :: this
    real(dp), intent(in) :: Emin, deltaE
    integer, intent(in) :: nE_global, nE_local
    real(dp), intent(in) :: kweights(:)
    integer, intent(in) :: local_kindex(:)

    integer :: iK, s, n_ev_out, n_ec_out
    real(dp) :: sumw, x, Emax

    if (deltaE <= 0.0_dp) call error_msg('SRH defect: energy step must be > 0')
    if (this%n_sites < 1) call error_msg('SRH defect: no defect site registered')

    this%Emin = Emin
    this%dE = deltaE
    this%nE = nE_global
    this%nE_local = nE_local
    this%iE_first = 1
    Emax = Emin + real(nE_global - 1, dp) * deltaE

    ! one rounded phonon shift used for EVERYTHING (ladder energies,
    ! Boltzmann factors, Bose factor, calibration)
    this%iEshift = max(1, nint(this%hw_in / deltaE))
    this%sh = real(this%iEshift, dp) * deltaE

    x = min(this%sh / this%kT, 700.0_dp)
    this%fB = 1.0_dp / (exp(x) - 1.0_dp)
    if (this%fB < 1.0e-300_dp) this%fB = 1.0e-300_dp

    this%width = max(this%width_in, 0.5_dp * deltaE)

    if (allocated(this%kw)) deallocate(this%kw)
    allocate(this%kw(size(local_kindex)))
    sumw = sum(kweights)
    if (sumw <= 0.0_dp) call error_msg('SRH defect: k weights sum to zero')
    do iK = 1, size(local_kindex)
      this%kw(iK) = kweights(local_kindex(iK)) / sumw
    end do

    call build_shapes(this)

    n_ev_out = 0
    n_ec_out = 0
    do s = 1, this%n_sites
      call setup_site(this, s)
      if (this%sites(s)%Ev < Emin) n_ev_out = n_ev_out + 1
      if (this%sites(s)%Ec > Emax) n_ec_out = n_ec_out + 1
    end do
    this%grid_ready = .true.

    if (id0) then
      write(*,'(A)') ' SRH defect: setup'
      write(*,'(A,I6)') '   number of defect sites=', this%n_sites
      write(*,'(A,F9.5,A,F9.5,A,I4,A)') '   phonon energy input=', this%hw_in, &
          & ' eV, used=', this%sh, ' eV (', this%iEshift, ' grid steps)'
      write(*,'(A,F9.6,A,F9.6,A)') '   energy step=', this%dE, ' eV, defect line width used=', &
          & this%width, ' eV'
      write(*,'(A,F10.5,A,F10.5,A)') '   Ed(abs) of first site=', this%sites(1)%Ed_abs, ' eV (Ev=', &
          & this%sites(1)%Ev, ' eV)'
      if (abs(this%sh - this%hw_in) > 0.1_dp * this%hw_in) then
        write(*,'(A)') '   WARNING: energy step too coarse: phonon energy rounded by more than 10%'
      end if
      if (this%width_in < 0.5_dp * deltaE) then
        write(*,'(A)') '   WARNING: energy step coarser than twice the defect width: accuracy of SRH rates is poor'
      end if
      if (this%width > 0.1_dp * this%kT) then
        write(*,'(A)') '   WARNING: defect line width > 0.1 kT: detailed balance is only approximate (error ~ width^2)'
      end if
      if (n_ev_out > 0) then
        write(*,'(A,I6,A)') '   WARNING: ', n_ev_out, ' site(s) have Ev below the energy window: hole capture is truncated'
      end if
      if (n_ec_out > 0) then
        write(*,'(A,I6,A)') '   WARNING: ', n_ec_out, ' site(s) have Ec above the energy window: electron capture is truncated'
      end if
    end if

  end subroutine set_EnGrid

#:if defined("MPI")
  !--------------------------------------------------------------------------
  !> Store the communicator topology and cache offset supplied by the NEGF
  !! driver.  The k communicator combines weighted local-k contributions at a
  !! fixed energy block; the energy communicator then distributes those blocks
  !! to every rank required by the phonon ladder.
  subroutine set_mpi_communicators(this, energy_comm, k_comm, iE_first)
    class(TSrhMultiphonon), intent(inout) :: this
    type(mpifx_comm), intent(in) :: energy_comm, k_comm
    integer, intent(in) :: iE_first

    this%energy_comm = energy_comm
    this%k_comm = k_comm
    this%iE_first = iE_first
  end subroutine set_mpi_communicators
#:endif

  !--------------------------------------------------------------------------
  !> Density of states prefactor of a parabolic band: N(E) = c sqrt(E_eV)
  !! [cm^-3 eV^-1], spin degeneracy 2 included.
  pure function dos_coef(m_eff) result(c)
    real(dp), intent(in) :: m_eff
    real(dp) :: c, c_si
    c_si = (1.0_dp / (2.0_dp * PI**2)) * (2.0_dp * m_eff * M0_KG / HBAR_J_S**2)**1.5_dp
    c = c_si * sqrt(Q_C) * Q_C * 1.0e-6_dp
  end function dos_coef

  !--------------------------------------------------------------------------
  !> log of the modified Bessel function I_l(x), series (x > 0)
  pure function log_bessel_i(l, x) result(lb)
    integer, intent(in) :: l
    real(dp), intent(in) :: x
    real(dp) :: lb
    integer, parameter :: MMAX = 300
    real(dp) :: t(0:MMAX-1), mx
    integer :: m
    do m = 0, MMAX - 1
      t(m) = real(2*m + l, dp) * log(0.5_dp * x) - log_gamma(real(m + 1, dp)) &
          & - log_gamma(real(m + l + 1, dp))
    end do
    mx = maxval(t)
    lb = mx + log(sum(exp(t - mx)))
  end function log_bessel_i

  !--------------------------------------------------------------------------
  !> Finite-temperature ladder shapes ([P] Eq. 6), common to all sites.
  !! Rescaled to max = 1; the overall constant is fixed per site by the
  !! calibration in setup_site.
  subroutine build_shapes(this)
    type(TSrhMultiphonon), intent(inout) :: this

    integer :: l, nl
    real(dp) :: x, lg, dn, dh, mxe, mxh
    real(dp), allocatable :: lse(:), lsh(:)

    nl = LMAX_RUNGS
    allocate(lse(0:nl), lsh(0:nl))
    if (allocated(this%shape_e)) deallocate(this%shape_e)
    if (allocated(this%shape_h)) deallocate(this%shape_h)
    if (allocated(this%bolt)) deallocate(this%bolt)
    allocate(this%shape_e(0:nl), this%shape_h(0:nl), this%bolt(0:nl))

    x = 2.0_dp * this%S_hr * sqrt(this%fB * (this%fB + 1.0_dp))
    do l = 0, nl
      lg = -this%S_hr * (2.0_dp * this%fB + 1.0_dp) &
          & + 0.5_dp * real(l, dp) * log((this%fB + 1.0_dp) / this%fB) &
          & + log_bessel_i(l, x)
      dn = (real(l, dp) - this%S_hr)**2 / this%S_hr
      dh = (real(l, dp) + this%S_hr)**2 / this%S_hr
      if (dn > 1.0e-300_dp) then
        lse(l) = log(dn) + lg
      else
        lse(l) = -1.0e300_dp
      end if
      lsh(l) = log(dh) + lg
      this%bolt(l) = exp(-real(l, dp) * this%sh / this%kT)
    end do
    mxe = maxval(lse)
    mxh = maxval(lsh)
    this%shape_e = exp(lse - mxe)
    this%shape_h = exp(lsh - mxh)

    deallocate(lse, lsh)

  end subroutine build_shapes

  !--------------------------------------------------------------------------
  !> Per-site setup: site check, coupling calibration ([A] App. B), defect line,
  !! storage of the band self-energies.
  subroutine setup_site(this, s)
    type(TSrhMultiphonon), intent(inout) :: this
    integer, intent(in) :: s

    integer :: l, nl, j, n, pl_start, pl_end, nrow
    real(dp) :: E, xx, dosc, dosv, Nc_eff, Nv_eff, den_n, den_p, kn, kp, pref, xj, sw
    character(len=200) :: msg

    nl = LMAX_RUNGS

    ! ---- site must lie inside its principal layer -------------------------
    if (this%sites(s)%pl > this%struct%num_PLs) then
      write(msg, '(A,I8,A,I8)') 'SRH defect: site ', s, ' has pl_index beyond the number of PLs ', this%struct%num_PLs
      call error_msg(trim(msg))
    end if
    pl_start = this%struct%mat_PL_start(this%sites(s)%pl)
    pl_end = this%struct%mat_PL_end(this%sites(s)%pl)
    nrow = pl_end - pl_start + 1
    if (this%sites(s)%orb + this%sites(s)%n_orb - 1 > nrow) then
      write(msg, '(A,I8,A)') 'SRH defect: site ', s, ' does not fit inside its principal layer (orb_index + n_orb)'
      call error_msg(trim(msg))
    end if

    ! ---- calibration with analytic parabolic bands, energies relative to Ev = 0:
    !   c_n = sum_l kappa(l) N_c(E_l) exp(-(E_l-Ec)/kT) / N_c,eff    ([A] Eq. B1)
    dosc = dos_coef(this%mass_c)
    dosv = dos_coef(this%mass_v)
    Nc_eff = dosc * 0.5_dp * sqrt(PI) * this%kT**1.5_dp
    Nv_eff = dosv * 0.5_dp * sqrt(PI) * this%kT**1.5_dp

    den_n = 0.0_dp
    do l = 0, nl
      E = this%Ed_rel + real(l, dp) * this%sh
      if (E <= this%sites(s)%Eg) cycle
      xx = E - this%sites(s)%Eg
      if (xx > 60.0_dp * this%kT) exit
      den_n = den_n + this%shape_e(l) * dosc * sqrt(xx) * exp(-xx / this%kT)
    end do
    den_p = 0.0_dp
    do l = 0, nl
      E = this%Ed_rel - real(l, dp) * this%sh
      if (E >= 0.0_dp) cycle
      xx = -E
      if (xx > 60.0_dp * this%kT) exit
      den_p = den_p + this%shape_h(l) * dosv * sqrt(xx) * exp(-xx / this%kT)
    end do
    if (den_n <= 0.0_dp .or. den_p <= 0.0_dp) then
      write(msg, '(A,I8,A)') 'SRH defect: site ', s, ': no phonon-ladder rung falls inside the conduction/valence ' // &
          & 'band within 60 kT; check phonon_frequency, defect_energy and band edges'
      call error_msg(trim(msg))
    end if

    kn = this%Cn * Nc_eff / den_n      ! [eV cm^3 / s]
    kp = this%Cp * Nv_eff / den_p

    ! M(l) = hbar * kappa(l) * g_s / (2 pi V_cell)  [eV^2]
    pref = HBAR_EV_S * G_SPIN / (2.0_dp * PI * this%cell_vol)
    if (allocated(this%sites(s)%Mn)) deallocate(this%sites(s)%Mn)
    if (allocated(this%sites(s)%Mp)) deallocate(this%sites(s)%Mp)
    allocate(this%sites(s)%Mn(0:nl), this%sites(s)%Mp(0:nl))
    do l = 0, nl
      this%sites(s)%Mn(l) = pref * kn * this%shape_e(l)
      this%sites(s)%Mp(l) = pref * kp * this%shape_h(l)
    end do

    ! ---- Gaussian representation of the delta-like defect spectral function;
    !      weights sum to 1
    xj = (this%sites(s)%Ed_abs - this%Emin) / this%dE + 1.0_dp
    if (xj < 1.0_dp .or. xj > real(this%nE, dp)) then
      write(msg, '(A,I8,A)') 'SRH defect: site ', s, ': defect energy is outside the NEGF energy window (Emin..Emax)'
      call error_msg(trim(msg))
    end if
    this%sites(s)%jlo = max(1, floor(xj - LINE_SIGMAS * this%width / this%dE))
    this%sites(s)%jhi = min(this%nE, ceiling(xj + LINE_SIGMAS * this%width / this%dE))
    n = this%sites(s)%jhi - this%sites(s)%jlo + 1
    if (n < 1) call error_msg('SRH defect: empty defect line')
    if (allocated(this%sites(s)%aline)) deallocate(this%sites(s)%aline)
    allocate(this%sites(s)%aline(n))
    sw = 0.0_dp
    do j = this%sites(s)%jlo, this%sites(s)%jhi
      E = this%Emin + real(j - 1, dp) * this%dE
      this%sites(s)%aline(j - this%sites(s)%jlo + 1) = exp(-0.5_dp * ((E - this%sites(s)%Ed_abs) / this%width)**2)
      sw = sw + this%sites(s)%aline(j - this%sites(s)%jlo + 1)
    end do
    if (sw <= 0.0_dp) call error_msg('SRH defect: defect line weights vanish')
    this%sites(s)%aline = this%sites(s)%aline / sw

    ! ---- storage and state --------------------------------------------------
    if (allocated(this%sites(s)%Sn_b)) deallocate(this%sites(s)%Sn_b)
    if (allocated(this%sites(s)%Sp_b)) deallocate(this%sites(s)%Sp_b)
    if (allocated(this%sites(s)%Sr_b)) deallocate(this%sites(s)%Sr_b)
    allocate(this%sites(s)%Sn_b(this%nE), this%sites(s)%Sp_b(this%nE), this%sites(s)%Sr_b(this%nE))
    if (allocated(this%sites(s)%Sn_bc)) deallocate(this%sites(s)%Sn_bc)
    if (allocated(this%sites(s)%Sp_bc)) deallocate(this%sites(s)%Sp_bc)
    allocate(this%sites(s)%Sn_bc(this%nE), this%sites(s)%Sp_bc(this%nE))
    this%sites(s)%Sn_bc = 0.0_dp
    this%sites(s)%Sp_bc = 0.0_dp
    this%sites(s)%rb_c = 0.0_dp
    this%sites(s)%rb_v = 0.0_dp
    this%sites(s)%res_band = 0.0_dp
    this%sites(s)%Sn_b = 0.0_dp
    this%sites(s)%Sp_b = 0.0_dp
    this%sites(s)%Sr_b = 0.0_dp
    this%sites(s)%fd = 0.5_dp
    this%sites(s)%omfd = 0.5_dp
    this%sites(s)%R_def = 0.0_dp
    this%sites(s)%U_rec = 0.0_dp

  end subroutine setup_site

  !--------------------------------------------------------------------------
  subroutine add_sigma_r_sp(this, esh, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(c_dns), dimension(:,:), intent(inout) :: esh
    integer, intent(in), optional :: en_index, k_index, spin
    ! single precision path is not used by TiberCAD (double precision caches)
  end subroutine add_sigma_r_sp

  !> esh = E S - H - Sigma^r, with Sigma^r = Sr_b - (i/2)(Sn_b + Sp_b)
  !! on the local site rows.
  subroutine add_sigma_r_dp(this, esh, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(z_dns), dimension(:,:), intent(inout) :: esh
    integer, intent(in), optional :: en_index, k_index, spin
    integer :: s, io, pl, iEg
    real(dp) :: gam, sr

    if (this%scba_iter .eq. 0) return
    if (.not. this%grid_ready) return
    if (.not. present(en_index)) return
    if (en_index < 1 .or. en_index > this%nE_local) return
    iEg = this%iE_first + en_index - 1
    if (iEg < 1 .or. iEg > this%nE) return

    do s = 1, this%n_sites
      gam = this%sites(s)%Sn_b(iEg) + this%sites(s)%Sp_b(iEg)
      sr = this%sites(s)%Sr_b(iEg)
      if (gam <= 0.0_dp .and. abs(sr) <= tiny(1.0_dp)) cycle
      pl = this%sites(s)%pl
      if (.not. allocated(esh(pl, pl)%val)) cycle
      do io = this%sites(s)%orb, this%sites(s)%orb + this%sites(s)%n_orb - 1
        esh(pl, pl)%val(io, io) = esh(pl, pl)%val(io, io) + cmplx(-sr, 0.5_dp * gam, dp)
      end do
    end do

  end subroutine add_sigma_r_dp

  subroutine add_sigma_n_sp(this, sigma, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(c_dns), dimension(:,:), intent(inout) :: sigma
    integer, intent(in), optional :: en_index, k_index, spin
    ! single precision path is not used by TiberCAD (double precision caches)
  end subroutine add_sigma_n_sp

  subroutine add_sigma_n_dp(this, sigma, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(z_dns), dimension(:,:), intent(inout) :: sigma
    integer, intent(in), optional :: en_index, k_index, spin
    integer :: s, io, pl, iEg

    if (this%scba_iter .eq. 0) return
    if (.not. this%grid_ready) return
    if (.not. present(en_index)) return
    if (en_index < 1 .or. en_index > this%nE_local) return
    iEg = this%iE_first + en_index - 1
    if (iEg < 1 .or. iEg > this%nE) return

    do s = 1, this%n_sites
      if (this%sites(s)%Sn_b(iEg) <= 0.0_dp) cycle
      pl = this%sites(s)%pl
      if (.not. allocated(sigma(pl, pl)%val)) cycle
      do io = this%sites(s)%orb, this%sites(s)%orb + this%sites(s)%n_orb - 1
        sigma(pl, pl)%val(io, io) = sigma(pl, pl)%val(io, io) + cmplx(this%sites(s)%Sn_b(iEg), 0.0_dp, dp)
      end do
    end do

  end subroutine add_sigma_n_dp

  subroutine get_sigma_n_blk(this, blk_sigma_n, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(z_dns), dimension(:,:), intent(inout) :: blk_sigma_n
    integer, intent(in), optional :: en_index, k_index, spin

    ! the library does not call this routine; kept consistent with add_sigma_n
    call this%add_sigma_n_dp(blk_sigma_n, en_index, k_index, spin)

  end subroutine get_sigma_n_blk

  subroutine get_sigma_n_mat(this, sigma_n, ii, jj, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(z_dns), intent(inout) :: sigma_n
    integer, intent(in) :: ii, jj
    integer, intent(in), optional :: en_index, k_index, spin
    integer :: s, io, iEg

    if (this%scba_iter .eq. 0) return
    if (.not. this%grid_ready) return
    if (ii /= jj) return
    if (.not. present(en_index)) return
    if (en_index < 1 .or. en_index > this%nE_local) return
    iEg = this%iE_first + en_index - 1
    if (iEg < 1 .or. iEg > this%nE) return
    if (.not. allocated(sigma_n%val)) return

    do s = 1, this%n_sites
      if (this%sites(s)%pl /= ii) cycle
      do io = this%sites(s)%orb, this%sites(s)%orb + this%sites(s)%n_orb - 1
        if (io <= size(sigma_n%val, 1)) &
            & sigma_n%val(io, io) = sigma_n%val(io, io) + cmplx(this%sites(s)%Sn_b(iEg), 0.0_dp, dp)
      end do
    end do

  end subroutine get_sigma_n_mat

  subroutine set_Gr(this, Gr, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(z_dns), dimension(:,:), intent(in) :: Gr
    integer, intent(in), optional :: en_index, k_index, spin
  end subroutine set_Gr

  subroutine set_Gn(this, Gn, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(z_dns), dimension(:,:), intent(in) :: Gn
    integer, intent(in), optional :: en_index, k_index, spin
  end subroutine set_Gn

  !--------------------------------------------------------------------------
  !> k-averaged local G^n and spectral function A at every site, for every
  !! energy index.  Returns .false. if the caches are not usable.
  function gather_local(this, Gn_all, A_all) result(ok)
    class(TSrhMultiphonon), intent(in) :: this
    real(dp), intent(out) :: Gn_all(:,:), A_all(:,:)
    logical :: ok

    integer :: nE, nK, iE, iEg, iK, s, io, cur_pl, pl, o1, o2
    real(dp) :: wk
    type(TMatLabel) :: label
    type(z_dns), pointer :: pGn, pGr

    ok = .false.
    Gn_all = 0.0_dp
    A_all = 0.0_dp
    if (.not. associated(this%G_n)) return
    if (.not. associated(this%G_r)) return

    select type(pGcache => this%G_n)
    type is (TMatrixCacheMem)
      if (.not. pGcache%isInitialized) return
      nE = size(pGcache%MatArray, 1)
      nK = size(pGcache%MatArray, 2)
    class default
      return
    end select

    if (nE /= this%nE_local) call error_msg('SRH defect: energy cache size differs from local energy grid')
    if (nK /= size(this%kw)) call error_msg('SRH defect: k-point cache size differs from k weights')

    label%spin = 1
    do iK = 1, nK
      label%kpoint = iK
      do iE = 1, nE
        iEg = this%iE_first + iE - 1
        if (iEg < 1 .or. iEg > this%nE) cycle
        label%energy_point = iE
        cur_pl = -1
        nullify(pGn)
        nullify(pGr)
        do s = 1, this%n_sites
          pl = this%sites(s)%pl
          if (pl /= cur_pl) then
            ! fetch the PL block once for all consecutive sites in the same PL
            label%row_block = pl
            label%col_block = pl
            nullify(pGn)
            nullify(pGr)
            if (this%G_n%is_cached(label)) call this%G_n%retrieve_pointer(pGn, label)
            if (this%G_r%is_cached(label)) call this%G_r%retrieve_pointer(pGr, label)
            cur_pl = pl
          end if
          o1 = this%sites(s)%orb
          o2 = o1 + this%sites(s)%n_orb - 1
          wk = this%kw(iK) / real(this%sites(s)%n_orb, dp)
          if (associated(pGn)) then
            if (allocated(pGn%val)) then
              if (o2 <= size(pGn%val, 1)) then
                do io = o1, o2
                  Gn_all(iEg, s) = Gn_all(iEg, s) + wk * real(pGn%val(io, io), dp)
                end do
              end if
            end if
          end if
          if (associated(pGr)) then
            if (allocated(pGr%val)) then
              if (o2 <= size(pGr%val, 1)) then
                do io = o1, o2
                  A_all(iEg, s) = A_all(iEg, s) - 2.0_dp * wk * aimag(pGr%val(io, io))
                end do
              end if
            end if
          end if
        end do
      end do
    end do
#:if defined("MPI")
    ! Each rank owns a local k subset and a contiguous energy block.  First
    ! finish the k average at fixed energy, then make every energy block
    ! available for the shifted multiphonon ladder.  This preserves the
    ! globally normalized weights assigned in set_EnGrid.
    call mpifx_allreduceip(this%k_comm, Gn_all, MPI_SUM)
    call mpifx_allreduceip(this%k_comm, A_all, MPI_SUM)
    call mpifx_allreduceip(this%energy_comm, Gn_all, MPI_SUM)
    call mpifx_allreduceip(this%energy_comm, A_all, MPI_SUM)
#:endif
    ok = .true.

  end function gather_local

  !--------------------------------------------------------------------------
  !> Causal real part of a retarded self-energy from its broadening:
  !!
  !!   Lambda(E) = P integral dE' Gamma(E') / [2 pi (E-E')].
  !!
  !! Gamma is represented as constant inside every uniform energy cell.  The
  !! analytical cell integral avoids the singular grid-point contribution and
  !! is considerably more accurate than omitting j=i from a point sum.  The
  !! finite simulation window is the same approximation already made for all
  !! energy-domain NEGF quantities; contributions outside it are unavailable.
  subroutine causal_realpart(gamma, lambda)
    real(dp), intent(in) :: gamma(:)
    real(dp), intent(out) :: lambda(:)

    integer :: i, j, n
    real(dp) :: xplus, xminus

    n = size(gamma)
    if (size(lambda) /= n) call error_msg('SRH defect: causal self-energy array size mismatch')
    lambda = 0.0_dp
    do i = 1, n
      do j = 1, n
        if (j == i .or. gamma(j) <= 0.0_dp) cycle
        xplus = real(i - j, dp) + 0.5_dp
        xminus = real(i - j, dp) - 0.5_dp
        lambda(i) = lambda(i) + gamma(j) * log(abs(xplus / xminus)) / (2.0_dp * PI)
      end do
    end do
  end subroutine causal_realpart

  !--------------------------------------------------------------------------
  !> Steps 1-5 of the module header for every site.  Called once per SCBA iteration.
  subroutine compute_Sigma_n(this, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    integer, intent(in), optional :: en_index, k_index, spin

    integer :: nE, ns, s, jj, l, ie_up, ie_dn
    real(dp), allocatable :: Gn_all(:,:), A_all(:,:), Gn_loc(:), Gp_loc(:)
    real(dp) :: aj, gam, xd, Rc, Rv, gross, bal, bal_max
    real(dp) :: sum_c, sum_v, res_max, rb_scale
    integer :: ie, n_res
    logical :: has_prev
    real(dp) :: fd_min, fd_max, U_min, U_max

    if (.not. this%grid_ready) return
    if (this%n_sites < 1) return

    nE = this%nE
    ns = this%n_sites
    allocate(Gn_all(nE, ns), A_all(nE, ns), Gn_loc(nE), Gp_loc(nE))
    if (.not. gather_local(this, Gn_all, A_all)) then
      deallocate(Gn_all, A_all, Gn_loc, Gp_loc)
      return
    end if

    bal_max = 0.0_dp
    res_max = 0.0_dp
    n_res = 0

    do s = 1, ns
      ! G^n >= 0 and G^p = A - G^n >= 0 (numerical noise clipped)
      Gn_loc = max(Gn_all(:, s), 0.0_dp)
      Gp_loc = max(A_all(:, s) - Gn_loc, 0.0_dp)

      ! ---- step 2: defect self-energies, averaged over the defect line ------
      this%sites(s)%Sn_dc = 0.0_dp
      this%sites(s)%Sp_dc = 0.0_dp
      this%sites(s)%Sn_dv = 0.0_dp
      this%sites(s)%Sp_dv = 0.0_dp
      do jj = this%sites(s)%jlo, this%sites(s)%jhi
        aj = this%sites(s)%aline(jj - this%sites(s)%jlo + 1)
        do l = 0, LMAX_RUNGS
          ie_up = jj + l * this%iEshift
          ie_dn = jj - l * this%iEshift
          if (ie_up > nE .and. ie_dn < 1) exit
          if (ie_up <= nE) then
            this%sites(s)%Sn_dc = this%sites(s)%Sn_dc + aj * this%sites(s)%Mn(l) * Gn_loc(ie_up)
            this%sites(s)%Sp_dc = this%sites(s)%Sp_dc + aj * this%sites(s)%Mn(l) * this%bolt(l) * Gp_loc(ie_up)
          end if
          if (ie_dn >= 1) then
            this%sites(s)%Sn_dv = this%sites(s)%Sn_dv + aj * this%sites(s)%Mp(l) * this%bolt(l) * Gn_loc(ie_dn)
            this%sites(s)%Sp_dv = this%sites(s)%Sp_dv + aj * this%sites(s)%Mp(l) * Gp_loc(ie_dn)
          end if
        end do
      end do

      ! ---- step 3: occupation, Eq. 59 ----------------------------------------
      gam = this%sites(s)%Sn_dc + this%sites(s)%Sp_dc + this%sites(s)%Sn_dv + this%sites(s)%Sp_dv
      if (gam > 0.0_dp) then
        this%sites(s)%fd = (this%sites(s)%Sn_dc + this%sites(s)%Sn_dv) / gam
        this%sites(s)%omfd = (this%sites(s)%Sp_dc + this%sites(s)%Sp_dv) / gam
        ! ---- step 4: net recombination per defect, closed form --------------
        this%sites(s)%R_def = (this%sites(s)%Sn_dc * this%sites(s)%Sp_dv &
            & - this%sites(s)%Sp_dc * this%sites(s)%Sn_dv) / gam / HBAR_EV_S
      else
        this%sites(s)%fd = 0.5_dp
        this%sites(s)%omfd = 0.5_dp
        this%sites(s)%R_def = 0.0_dp
      end if
      this%sites(s)%U_rec = this%rho_d * this%sites(s)%R_def

      ! consistency diagnostic R_c = R_v (Eqs. 55, 58)
      Rc = (this%sites(s)%Sn_dc * this%sites(s)%omfd - this%sites(s)%Sp_dc * this%sites(s)%fd) / HBAR_EV_S
      Rv = (this%sites(s)%Sp_dv * this%sites(s)%fd - this%sites(s)%Sn_dv * this%sites(s)%omfd) / HBAR_EV_S
      gross = max(abs(this%sites(s)%Sn_dc * this%sites(s)%omfd), abs(this%sites(s)%Sp_dc * this%sites(s)%fd), &
          & abs(this%sites(s)%Sp_dv * this%sites(s)%fd), abs(this%sites(s)%Sn_dv * this%sites(s)%omfd)) / HBAR_EV_S
      if (gross > 0.0_dp) then
        bal = abs(Rc - Rv) / gross
      else
        bal = 0.0_dp
      end if
      bal_max = max(bal_max, bal)

      ! ---- step 5: band self-energies from the defect Green's function ------
      xd = this%rho_d * this%cell_vol * this%sites(s)%weight      ! defects per orbital cell

      ! ---- SCBA / conservation check: OLD band self-energies x NEW G -----------
      ! Net electron loss of the bands per defect, (1/(2 pi hbar x_d)) int dE [Sp Gn - Sn Gp],
      ! evaluated with the band self-energies stored from the previous iteration and
      ! the G just returned by the NEGF solve (which already contains them).
      !   rb_c : conduction band -> defect, rb_v : defect -> valence band.
      ! At the SCBA fixed point they equal the defect-side rates Rc, Rv computed above
      ! from the same G.  Their difference is therefore a real convergence/consistency
      ! measure: it tests that the band self-energies (2 pi, x_d, ladder directions,
      ! Boltzmann factors) really enter G and reproduce the defect balance.
      has_prev = any(this%sites(s)%Sn_b > 0.0_dp .or. this%sites(s)%Sp_b > 0.0_dp)
      this%sites(s)%rb_c = 0.0_dp
      this%sites(s)%rb_v = 0.0_dp
      if (has_prev .and. xd > 0.0_dp) then
        sum_c = 0.0_dp
        sum_v = 0.0_dp
        do ie = 1, nE
          sum_c = sum_c + this%sites(s)%Sp_bc(ie) * Gn_loc(ie) - this%sites(s)%Sn_bc(ie) * Gp_loc(ie)
          sum_v = sum_v + (this%sites(s)%Sn_b(ie) - this%sites(s)%Sn_bc(ie)) * Gp_loc(ie) &
              & - (this%sites(s)%Sp_b(ie) - this%sites(s)%Sp_bc(ie)) * Gn_loc(ie)
        end do
        rb_scale = this%dE / (2.0_dp * PI * HBAR_EV_S * xd)
        this%sites(s)%rb_c = sum_c * rb_scale
        this%sites(s)%rb_v = sum_v * rb_scale
      else
        has_prev = .false.
      end if

      this%sites(s)%Sn_b = 0.0_dp
      this%sites(s)%Sp_b = 0.0_dp
      this%sites(s)%Sn_bc = 0.0_dp
      this%sites(s)%Sp_bc = 0.0_dp
      do jj = this%sites(s)%jlo, this%sites(s)%jhi
        ! A_d(E_j) = i(G^R_d - G^A_d) = 2 pi delta(E-E_d) ([A] Eq. 7, [P] Eq. 1-2): the
        ! spectral function integrates to 2 pi, not 1.  The defect-side sums (step 2) use
        ! the normalised line (sum = 1) because the 2 pi cancels the 1/(2 pi hbar) of
        ! Eq. 55/58; here the 2 pi must be explicit, otherwise Sigma_b is 2 pi too small
        ! and band/defect rates do not balance.
        aj = 2.0_dp * PI * this%sites(s)%aline(jj - this%sites(s)%jlo + 1) / this%dE      ! A_d(E_j) [1/eV]
        do l = 0, LMAX_RUNGS
          ie_up = jj + l * this%iEshift      ! E = E_j + l hw : conduction-type
          ie_dn = jj - l * this%iEshift      ! E = E_j - l hw : valence-type
          if (ie_up > nE .and. ie_dn < 1) exit
          if (ie_up <= nE) then
            this%sites(s)%Sn_b(ie_up) = this%sites(s)%Sn_b(ie_up) &
                & + xd * this%sites(s)%Mn(l) * this%bolt(l) * this%sites(s)%fd * aj
            this%sites(s)%Sp_b(ie_up) = this%sites(s)%Sp_b(ie_up) &
                & + xd * this%sites(s)%Mn(l) * this%sites(s)%omfd * aj
            this%sites(s)%Sn_bc(ie_up) = this%sites(s)%Sn_bc(ie_up) &
                & + xd * this%sites(s)%Mn(l) * this%bolt(l) * this%sites(s)%fd * aj
            this%sites(s)%Sp_bc(ie_up) = this%sites(s)%Sp_bc(ie_up) &
                & + xd * this%sites(s)%Mn(l) * this%sites(s)%omfd * aj
          end if
          if (ie_dn >= 1) then
            this%sites(s)%Sn_b(ie_dn) = this%sites(s)%Sn_b(ie_dn) &
                & + xd * this%sites(s)%Mp(l) * this%sites(s)%fd * aj
            this%sites(s)%Sp_b(ie_dn) = this%sites(s)%Sp_b(ie_dn) &
                & + xd * this%sites(s)%Mp(l) * this%bolt(l) * this%sites(s)%omfd * aj
          end if
        end do
      end do
      ! Retarded and lesser/greater components must describe the same causal
      ! interaction.  Recover the principal-value part from Gamma=Sn+Sp.
      call causal_realpart(this%sites(s)%Sn_b + this%sites(s)%Sp_b, this%sites(s)%Sr_b)

      if (has_prev) then
        if (gross > 0.0_dp) then
          this%sites(s)%res_band = max(abs(this%sites(s)%rb_c - Rc), abs(this%sites(s)%rb_v - Rv)) / gross
        else
          this%sites(s)%res_band = 0.0_dp
        end if
        res_max = max(res_max, this%sites(s)%res_band)
        n_res = n_res + 1
      end if
    end do

    if (id0) then
      if (ns <= 3) then
        do s = 1, ns
          write(*,'(A,I3,A,I5,A,F9.6,A,ES11.3,A,ES11.3)') ' SRH defect (scba_iter=', this%scba_iter, &
              & ', site ', s, '): fd=', this%sites(s)%fd, '  R_def[1/s]=', this%sites(s)%R_def, &
              & '  U[cm-3 s-1]=', this%sites(s)%U_rec
        end do
      else
        fd_min = minval(this%sites(1:ns)%fd)
        fd_max = maxval(this%sites(1:ns)%fd)
        U_min = minval(this%sites(1:ns)%U_rec)
        U_max = maxval(this%sites(1:ns)%U_rec)
        write(*,'(A,I3,A,I6,A,F9.6,A,F9.6,A,ES10.3,A,ES10.3)') ' SRH defects (scba_iter=', this%scba_iter, &
            & ', ', ns, ' sites): fd in [', fd_min, ',', fd_max, ']  U[cm-3 s-1] in [', U_min, ',', U_max
      end if
      write(*,'(A,ES9.2)') '   max |Rc-Rv|/gross over sites = ', bal_max
      if (n_res > 0) then
        write(*,'(A,ES9.2)') '   SCBA band/defect flux residual |rb(old Sigma_b)-R(new G)|/gross, max over sites = ', res_max
      else
        write(*,'(A)') '   SCBA band/defect flux residual: n/a (first iteration, no previous Sigma_b)'
      end if
    end if

    deallocate(Gn_all, A_all, Gn_loc, Gp_loc)

  end subroutine compute_Sigma_n

  !--------------------------------------------------------------------------
  !> Sigma^r is built on the fly in add_sigma_r_dp from the causal real part
  !! and anti-Hermitian broadening prepared with the same lesser/greater
  !! self-energies in compute_Sigma_n.
  subroutine compute_Sigma_r(this, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    integer, intent(in), optional :: en_index, k_index, spin
  end subroutine compute_Sigma_r

  !--------------------------------------------------------------------------
  !> The self-energies are overwritten by compute_Sigma_n every iteration; they
  !! must survive until then (compute_Sigma_r / add_sigma_* read them), so
  !! nothing is deleted here.
  subroutine destroy_Sigma_r(this)
    class(TSrhMultiphonon) :: this
  end subroutine destroy_Sigma_r

  subroutine destroy_Sigma_n(this)
    class(TSrhMultiphonon) :: this
  end subroutine destroy_Sigma_n

  subroutine destroy(this)
    class(TSrhMultiphonon) :: this
    if (associated(this%sigma_r)) then
      call this%sigma_r%destroy()
      deallocate(this%sigma_r)
      nullify(this%sigma_r)
    end if
    if (associated(this%sigma_n)) then
      call this%sigma_n%destroy()
      deallocate(this%sigma_n)
      nullify(this%sigma_n)
    end if
    if (allocated(this%sites)) deallocate(this%sites)
    this%n_sites = 0
    if (allocated(this%shape_e)) deallocate(this%shape_e)
    if (allocated(this%shape_h)) deallocate(this%shape_h)
    if (allocated(this%bolt)) deallocate(this%bolt)
    if (allocated(this%kw)) deallocate(this%kw)
  end subroutine destroy

end module srh_multiphonon
