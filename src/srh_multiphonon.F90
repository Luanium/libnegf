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
!! One localized defect level (delta-like spectral function, Eq. 8 of [A]) !
!! at the site (pl_index; orbitals orb_index .. orb_index+n_orb-1).        !
!! Per SCBA iteration:                                                      !
!!  1. local band Green's functions G^n, G^p = A - G^n at the defect site    !
!!     (average over the n_orb orbitals of the site, then over k)           !
!!     (k-averaged) are read from the libNEGF caches (full G: conduction   !
!!     and valence band are both contained, no separation needed).          !
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
!!     G^n_d = fd A_d, G^p_d = (1-fd) A_d, written on EVERY orbital of the  !
!!     defect site (same value; flux-consistent, see Limitations):         !
!!       S^n_cd(E) = x_d sum_l Mn(l) e^{-b l hw} fd     A_d(E - l hw)       !
!!       S^p_cd(E) = x_d sum_l Mn(l)             (1-fd) A_d(E - l hw)       !
!!       S^n_vd(E) = x_d sum_l Mp(l)             fd     A_d(E + l hw)       !
!!       S^p_vd(E) = x_d sum_l Mp(l) e^{-b l hw} (1-fd) A_d(E + l hw)       !
!!     S^r = -(i/2)(S^n + S^p)   (anti-Hermitian part, level shift dropped) !
!!                                                                          !
!! Conventions: libNEGF G_n = -i G^< (Hermitian, real positive diagonal),   !
!! G_p = i G^> = A - G_n, A = -2 Im G_r.  Sigma_n, Sigma_p likewise.        !
!!                                                                          !
!! Coupling strength: M(l) = hbar kappa shape(l) g_s / (2 pi V_cell), with  !
!! shape(l) the finite-temperature ladder weight of [P] Eq. 6 and kappa     !
!! fixed by the thermal-average calibration of [A] App. B (Eq. B1/B2) with  !
!! analytic parabolic bands: c_n = Cn = sigma_n v_th (c_p = Cp likewise).   !
!! Hole weight uses (l+S)^2/S, electron weight (l-S)^2/S ([A] Eq. A5, A10).  !
!!                                                                          !
!! Limitations (first complete version): one MPI rank on energy and k axes; !
!! one defect per SITE: in a multi-orbital basis the conduction and valence !
!! states sit on different orbitals, so the defect sees the orbital average !
!! G_loc = (1/n_orb) sum_o G_oo and the same Sigma is put on every orbital: !
!! the exchanged flux is then sum_o Sigma G_oo = n_orb Sigma G_loc, which   !
!! equals rho_d V_site R_def when cell_volume is the volume per ORBITAL     !
!! (V_site / n_orb).  For n_orb = 1 this is the single-orbital model.       !
!! Uniform energy grid; defect energy fixed in                              !
!! absolute energy (local band edges given by the user).                    !
!!--------------------------------------------------------------------------!

module srh_multiphonon

  use ln_precision, only : dp
  use ln_messages, only : error_msg
  use mpi_globals, only : id0
  use interactions, only : TInteraction
  use ln_inelastic, only : TInelastic
  use ln_structure, only : TStruct_info
  use mat_def, only : c_dns, z_dns, create, mat_destroy => destroy
  use ln_cache, only : TMatrixCacheMem, TMatLabel

  implicit none
  private

  public :: TSrhMultiphonon
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

  type, extends(TInelastic) :: TSrhMultiphonon

    ! ---- user parameters -------------------------------------------------
    !> defect level above the local valence band edge [eV]
    real(dp) :: Ed_rel = 0.0_dp
    !> absolute energy of the defect level on the NEGF energy axis [eV]
    real(dp) :: Ed_abs = 0.0_dp
    !> band edges at the defect site, absolute energies [eV]
    real(dp) :: Ev_abs = 0.0_dp
    real(dp) :: Ec_abs = 0.0_dp
    !> gap [eV]
    real(dp) :: Eg = 0.0_dp
    !> defect density [cm^-3]
    real(dp) :: rho_d = 0.0_dp
    !> capture coefficients sigma*v_th [cm^3/s]
    real(dp) :: Cn = 0.0_dp
    real(dp) :: Cp = 0.0_dp
    real(dp) :: temperature = 300.0_dp
    !> defect site: principal layer, first orbital (row inside the PL block, 1-based)
    !! and number of consecutive orbitals that belong to the site
    integer :: pl_index = 0
    integer :: orb_index = 0
    integer :: n_orb = 1
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

    ! ---- derived: energy grid -------------------------------------------
    logical :: grid_ready = .false.
    real(dp) :: Emin = 0.0_dp
    real(dp) :: dE = 0.0_dp
    integer :: nE = 0
    integer :: iEshift = 0
    !> rounded phonon energy actually used everywhere (iEshift*dE) [eV]
    real(dp) :: sh = 0.0_dp
    real(dp) :: kT = 0.0_dp
    real(dp) :: fB = 0.0_dp
    real(dp) :: width = 0.0_dp
    real(dp), allocatable :: kw(:)

    ! ---- derived: ladder -------------------------------------------------
    !> M_dc(l), M_dv(l) [eV^2]
    real(dp), allocatable :: Mn(:)
    real(dp), allocatable :: Mp(:)
    !> Boltzmann factor exp(-l sh / kT)
    real(dp), allocatable :: bolt(:)

    ! ---- derived: defect line --------------------------------------------
    integer :: jlo = 0
    integer :: jhi = 0
    real(dp), allocatable :: aline(:)

    ! ---- state -----------------------------------------------------------
    real(dp) :: fd = 0.5_dp
    real(dp) :: omfd = 0.5_dp
    !> net recombination per defect [1/s] and volume rate U = rho_d R [cm^-3 s^-1]
    real(dp) :: R_def = 0.0_dp
    real(dp) :: U_rec = 0.0_dp
    !> line-averaged defect self-energies [eV]
    real(dp) :: Sn_dc = 0.0_dp
    real(dp) :: Sp_dc = 0.0_dp
    real(dp) :: Sn_dv = 0.0_dp
    real(dp) :: Sp_dv = 0.0_dp
    !> band self-energies at the defect orbital as function of energy index [eV]
    real(dp), allocatable :: Sn_b(:)
    real(dp), allocatable :: Sp_b(:)

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
    procedure :: set_EnGrid

  end type TSrhMultiphonon

contains

  !--------------------------------------------------------------------------
  subroutine TSrhMultiphonon_create(this)
    class(TInteraction), allocatable :: this
    allocate(TSrhMultiphonon::this)
  end subroutine TSrhMultiphonon_create

  !--------------------------------------------------------------------------
  !> Store the physical parameters (no energy grid yet; see set_EnGrid)
  subroutine TSrhMultiphonon_init(this, struct, Ed_rel, rho_d, Cn, Cp, temperature, &
      & pl_index, orb_index, niter, phonon_energy, huang_rhys, mass_c, mass_v, &
      & Ev_abs, Ec_abs, cell_volume, line_width)
    type(TSrhMultiphonon), intent(inout) :: this
    type(TStruct_info), intent(in) :: struct
    real(dp), intent(in) :: Ed_rel, rho_d, Cn, Cp, temperature
    integer, intent(in) :: pl_index, orb_index, niter
    real(dp), intent(in) :: phonon_energy, huang_rhys, mass_c, mass_v
    real(dp), intent(in) :: Ev_abs, Ec_abs, cell_volume, line_width

    if (phonon_energy <= 0.0_dp) call error_msg('SRH defect: phonon_frequency must be > 0 (eV)')
    if (huang_rhys <= 0.0_dp) call error_msg('SRH defect: huang_rhys must be > 0')
    if (mass_c <= 0.0_dp .or. mass_v <= 0.0_dp) &
        & call error_msg('SRH defect: effective_mass_c and effective_mass_v must be > 0 (units of m0)')
    if (Ec_abs - Ev_abs <= 0.0_dp) &
        & call error_msg('SRH defect: band_edge_c must be above band_edge_v (eV, absolute)')
    if (Ed_rel <= 0.0_dp .or. Ed_rel >= Ec_abs - Ev_abs) &
        & call error_msg('SRH defect: defect_energy (above Ev) must lie inside the gap')
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
    this%pl_index = pl_index
    this%orb_index = orb_index
    this%scba_niter = niter
    this%hw_in = phonon_energy
    this%S_hr = huang_rhys
    this%mass_c = mass_c
    this%mass_v = mass_v
    this%Ev_abs = Ev_abs
    this%Ec_abs = Ec_abs
    this%Eg = Ec_abs - Ev_abs
    this%Ed_abs = Ev_abs + Ed_rel
    this%cell_vol = cell_volume
    if (line_width > 0.0_dp) this%width_in = line_width
    this%kT = KB_EV * temperature

    this%fd = 0.5_dp
    this%omfd = 0.5_dp
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
  !> Called by libNEGF before the SCBA loop. Builds everything that depends on
  !! the energy grid: rounded phonon shift, ladder weights, calibration,
  !! defect line, normalised k weights.
  subroutine set_EnGrid(this, Emin, deltaE, nE_global, nE_local, kweights, local_kindex)
    class(TSrhMultiphonon), intent(inout) :: this
    real(dp), intent(in) :: Emin, deltaE
    integer, intent(in) :: nE_global, nE_local
    real(dp), intent(in) :: kweights(:)
    integer, intent(in) :: local_kindex(:)

    integer :: iK
    real(dp) :: sumw, x

    if (deltaE <= 0.0_dp) call error_msg('SRH defect: energy step must be > 0')
    if (nE_local /= nE_global) &
        & call error_msg('SRH defect: parallelization over the energy axis is not supported (use 1 energy rank)')
    if (size(local_kindex) /= size(kweights)) &
        & call error_msg('SRH defect: parallelization over k points is not supported (use 1 k rank)')

    this%Emin = Emin
    this%dE = deltaE
    this%nE = nE_global

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

    call build_ladder_and_calibrate(this)
    call build_line(this)

    if (allocated(this%Sn_b)) deallocate(this%Sn_b)
    if (allocated(this%Sp_b)) deallocate(this%Sp_b)
    allocate(this%Sn_b(this%nE), this%Sp_b(this%nE))
    this%Sn_b = 0.0_dp
    this%Sp_b = 0.0_dp
    this%grid_ready = .true.

    if (id0) then
      write(*,'(A)') ' SRH defect: setup'
      write(*,'(A,F9.5,A,F9.5,A,I4,A)') '   phonon energy input=', this%hw_in, &
          & ' eV, used=', this%sh, ' eV (', this%iEshift, ' grid steps)'
      write(*,'(A,F9.6,A,F9.6,A,I4,A)') '   energy step=', this%dE, ' eV, defect line width used=', &
          & this%width, ' eV (', this%jhi - this%jlo + 1, ' grid points)'
      write(*,'(A,F10.5,A,F10.5,A,F10.5)') '   Ed(abs)=', this%Ed_abs, ' eV, Ev=', this%Ev_abs, ' eV, Ec=', this%Ec_abs
      if (abs(this%sh - this%hw_in) > 0.1_dp * this%hw_in) then
        write(*,'(A)') '   WARNING: energy step too coarse: phonon energy rounded by more than 10%'
      end if
      if (this%width_in < 0.5_dp * deltaE) then
        write(*,'(A)') '   WARNING: energy step coarser than twice the defect width: accuracy of SRH rates is poor'
      end if
      if (this%width > 0.1_dp * this%kT) then
        write(*,'(A)') '   WARNING: defect line width > 0.1 kT: detailed balance is only approximate (error ~ width^2)'
      end if
    end if

  end subroutine set_EnGrid

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
  !> Ladder weights ([P] Eq. 6) and calibration of the coupling ([A] App. B).
  subroutine build_ladder_and_calibrate(this)
    type(TSrhMultiphonon), intent(inout) :: this

    integer :: l, nl
    real(dp) :: x, lg, dn, dh, mxe, mxh, E, xx
    real(dp) :: dosc, dosv, Nc_eff, Nv_eff, den_n, den_p, kn, kp, pref
    real(dp), allocatable :: lse(:), lsh(:), wsh_e(:), wsh_h(:)

    nl = LMAX_RUNGS
    allocate(lse(0:nl), lsh(0:nl), wsh_e(0:nl), wsh_h(0:nl))
    if (allocated(this%Mn)) deallocate(this%Mn)
    if (allocated(this%Mp)) deallocate(this%Mp)
    if (allocated(this%bolt)) deallocate(this%bolt)
    allocate(this%Mn(0:nl), this%Mp(0:nl), this%bolt(0:nl))

    ! finite-temperature ladder shape (log domain, rescaled to max = 1;
    ! the overall constant is fixed by the calibration below)
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
    wsh_e = exp(lse - mxe)
    wsh_h = exp(lsh - mxh)

    ! calibration with analytic parabolic bands, energies relative to Ev = 0:
    !   c_n = sum_l kappa(l) N_c(E_l) exp(-(E_l-Ec)/kT) / N_c,eff    ([A] Eq. B1)
    dosc = dos_coef(this%mass_c)
    dosv = dos_coef(this%mass_v)
    Nc_eff = dosc * 0.5_dp * sqrt(PI) * this%kT**1.5_dp
    Nv_eff = dosv * 0.5_dp * sqrt(PI) * this%kT**1.5_dp

    den_n = 0.0_dp
    do l = 0, nl
      E = this%Ed_rel + real(l, dp) * this%sh
      if (E <= this%Eg) cycle
      xx = E - this%Eg
      if (xx > 60.0_dp * this%kT) exit
      den_n = den_n + wsh_e(l) * dosc * sqrt(xx) * exp(-xx / this%kT)
    end do
    den_p = 0.0_dp
    do l = 0, nl
      E = this%Ed_rel - real(l, dp) * this%sh
      if (E >= 0.0_dp) cycle
      xx = -E
      if (xx > 60.0_dp * this%kT) exit
      den_p = den_p + wsh_h(l) * dosv * sqrt(xx) * exp(-xx / this%kT)
    end do
    if (den_n <= 0.0_dp .or. den_p <= 0.0_dp) then
      call error_msg('SRH defect: no phonon-ladder rung falls inside the conduction/valence band ' // &
          & 'within 60 kT; check phonon_frequency, defect_energy and band edges')
    end if

    kn = this%Cn * Nc_eff / den_n      ! [eV cm^3 / s]
    kp = this%Cp * Nv_eff / den_p

    ! M(l) = hbar * kappa(l) * g_s / (2 pi V_cell)  [eV^2]
    pref = HBAR_EV_S * G_SPIN / (2.0_dp * PI * this%cell_vol)
    do l = 0, nl
      this%Mn(l) = pref * kn * wsh_e(l)
      this%Mp(l) = pref * kp * wsh_h(l)
    end do

    deallocate(lse, lsh, wsh_e, wsh_h)

  end subroutine build_ladder_and_calibrate

  !--------------------------------------------------------------------------
  !> Gaussian representation of the delta-like defect spectral function on the
  !! energy grid; weights sum to 1.
  subroutine build_line(this)
    type(TSrhMultiphonon), intent(inout) :: this

    integer :: j, n
    real(dp) :: xj, E, s

    xj = (this%Ed_abs - this%Emin) / this%dE + 1.0_dp
    if (xj < 1.0_dp .or. xj > real(this%nE, dp)) then
      call error_msg('SRH defect: defect energy is outside the NEGF energy window (Emin..Emax)')
    end if
    this%jlo = max(1, floor(xj - LINE_SIGMAS * this%width / this%dE))
    this%jhi = min(this%nE, ceiling(xj + LINE_SIGMAS * this%width / this%dE))
    n = this%jhi - this%jlo + 1
    if (n < 1) call error_msg('SRH defect: empty defect line')
    if (allocated(this%aline)) deallocate(this%aline)
    allocate(this%aline(n))
    s = 0.0_dp
    do j = this%jlo, this%jhi
      E = this%Emin + real(j - 1, dp) * this%dE
      this%aline(j - this%jlo + 1) = exp(-0.5_dp * ((E - this%Ed_abs) / this%width)**2)
      s = s + this%aline(j - this%jlo + 1)
    end do
    if (s <= 0.0_dp) call error_msg('SRH defect: defect line weights vanish')
    this%aline = this%aline / s

  end subroutine build_line

  !--------------------------------------------------------------------------
  logical function site_ok(this)
    class(TSrhMultiphonon), intent(in) :: this
    integer :: nbl, pl_start, pl_end, nrow

    site_ok = .false.
    if (this%pl_index < 1 .or. this%orb_index < 1 .or. this%n_orb < 1) return
    nbl = this%struct%num_PLs
    if (this%pl_index > nbl) return
    pl_start = this%struct%mat_PL_start(this%pl_index)
    pl_end   = this%struct%mat_PL_end(this%pl_index)
    nrow = pl_end - pl_start + 1
    ! the whole site (n_orb consecutive rows) must lie inside the PL block
    if (this%orb_index + this%n_orb - 1 > nrow) return
    site_ok = .true.
  end function site_ok

  !--------------------------------------------------------------------------
  subroutine add_sigma_r_sp(this, esh, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(c_dns), dimension(:,:), intent(inout) :: esh
    integer, intent(in), optional :: en_index, k_index, spin
    ! single precision path is not used by TiberCAD (double precision caches)
  end subroutine add_sigma_r_sp

  subroutine add_sigma_r_dp(this, esh, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(z_dns), dimension(:,:), intent(inout) :: esh
    integer, intent(in), optional :: en_index, k_index, spin
    type(TMatLabel) :: label
    type(z_dns), pointer :: tmp
    integer :: pl

    if (this%scba_iter .eq. 0) return
    if (.not. site_ok(this)) return
    if (.not. associated(this%sigma_r)) return
    if (.not. present(en_index) .or. .not. present(k_index)) return

    pl = this%pl_index
    label%energy_point = en_index
    label%kpoint = k_index
    label%spin = 1
    label%row_block = pl
    label%col_block = pl

    if (.not. this%sigma_r%is_cached(label)) return
    call this%sigma_r%retrieve_pointer(tmp, label)
    esh(pl, pl)%val = esh(pl, pl)%val - tmp%val

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
    type(TMatLabel) :: label
    type(z_dns), pointer :: tmp
    integer :: pl

    if (this%scba_iter .eq. 0) return
    if (.not. site_ok(this)) return
    if (.not. associated(this%sigma_n)) return
    if (.not. present(en_index) .or. .not. present(k_index)) return

    pl = this%pl_index
    label%energy_point = en_index
    label%kpoint = k_index
    label%spin = 1
    label%row_block = pl
    label%col_block = pl

    if (.not. this%sigma_n%is_cached(label)) return
    call this%sigma_n%retrieve_pointer(tmp, label)
    sigma(pl, pl)%val = sigma(pl, pl)%val + tmp%val

  end subroutine add_sigma_n_dp

  subroutine get_sigma_n_blk(this, blk_sigma_n, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(z_dns), dimension(:,:), intent(inout) :: blk_sigma_n
    integer, intent(in), optional :: en_index, k_index, spin
    type(TMatLabel) :: label
    integer :: pl

    if (this%scba_iter .eq. 0) return
    if (.not. site_ok(this)) return
    if (.not. associated(this%sigma_n)) return
    if (.not. present(en_index) .or. .not. present(k_index)) return

    pl = this%pl_index
    label%energy_point = en_index
    label%kpoint = k_index
    label%spin = 1
    label%row_block = pl
    label%col_block = pl
    if (this%sigma_n%is_cached(label)) then
      call this%sigma_n%retrieve(blk_sigma_n(pl, pl), label)
    end if

  end subroutine get_sigma_n_blk

  subroutine get_sigma_n_mat(this, sigma_n, ii, jj, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(z_dns), intent(inout) :: sigma_n
    integer, intent(in) :: ii, jj
    integer, intent(in), optional :: en_index, k_index, spin
    type(TMatLabel) :: label

    if (this%scba_iter .eq. 0) return
    if (.not. site_ok(this)) return
    if (ii /= this%pl_index .or. jj /= this%pl_index) return
    if (.not. associated(this%sigma_n)) return
    if (.not. present(en_index) .or. .not. present(k_index)) return

    label%energy_point = en_index
    label%kpoint = k_index
    label%spin = 1
    label%row_block = ii
    label%col_block = jj
    if (this%sigma_n%is_cached(label)) then
      call this%sigma_n%retrieve(sigma_n, label)
    end if

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
  !> k-averaged local G^n and spectral function A at the defect orbital,
  !! for every energy index.  Returns .false. if the caches are not usable.
  function gather_local(this, Gn_loc, A_loc) result(ok)
    class(TSrhMultiphonon), intent(in) :: this
    real(dp), intent(out) :: Gn_loc(:), A_loc(:)
    logical :: ok

    integer :: nE, nK, iE, iK, pl, orb
    type(TMatLabel) :: label
    type(z_dns), pointer :: pG

    ok = .false.
    Gn_loc = 0.0_dp
    A_loc = 0.0_dp
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

    if (nE /= this%nE) call error_msg('SRH defect: energy cache size differs from energy grid size')
    if (nK /= size(this%kw)) call error_msg('SRH defect: k-point cache size differs from k weights')

    pl = this%pl_index
    orb = this%orb_index
    label%spin = 1
    label%row_block = pl
    label%col_block = pl

    do iK = 1, nK
      do iE = 1, nE
        label%kpoint = iK
        label%energy_point = iE
        if (this%G_n%is_cached(label)) then
          call this%G_n%retrieve_pointer(pG, label)
          if (associated(pG)) then
            if (allocated(pG%val)) then
              if (orb + this%n_orb - 1 <= size(pG%val, 1)) then
                block
                  integer :: io
                  ! orbital average over the defect site, then k weight
                  do io = orb, orb + this%n_orb - 1
                    Gn_loc(iE) = Gn_loc(iE) + this%kw(iK) * real(pG%val(io, io), dp) &
                        & / real(this%n_orb, dp)
                  end do
                end block
              end if
            end if
          end if
        end if
        if (this%G_r%is_cached(label)) then
          call this%G_r%retrieve_pointer(pG, label)
          if (associated(pG)) then
            if (allocated(pG%val)) then
              if (orb + this%n_orb - 1 <= size(pG%val, 1)) then
                block
                  integer :: io
                  do io = orb, orb + this%n_orb - 1
                    A_loc(iE) = A_loc(iE) - 2.0_dp * this%kw(iK) * aimag(pG%val(io, io)) &
                        & / real(this%n_orb, dp)
                  end do
                end block
              end if
            end if
          end if
        end if
      end do
    end do
    ok = .true.

  end function gather_local

  !--------------------------------------------------------------------------
  !> Steps 1-5 of the module header.  Called once per SCBA iteration.
  subroutine compute_Sigma_n(this, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    integer, intent(in), optional :: en_index, k_index, spin

    integer :: nE, nK, nbl, pl, orb, nrow, iE, iK, jj, l, ie_up, ie_dn
    type(TMatLabel) :: label
    type(z_dns) :: sig
    real(dp), allocatable :: Gn_loc(:), Gp_loc(:), A_loc(:)
    real(dp) :: aj, gam, xd, Rc, Rv, gross, bal

    if (.not. site_ok(this)) return
    if (.not. this%grid_ready) return
    if (.not. associated(this%sigma_n)) return

    nE = this%nE
    allocate(Gn_loc(nE), Gp_loc(nE), A_loc(nE))
    if (.not. gather_local(this, Gn_loc, A_loc)) then
      deallocate(Gn_loc, Gp_loc, A_loc)
      return
    end if
    ! G^n >= 0 and G^p = A - G^n >= 0 (numerical noise clipped)
    Gn_loc = max(Gn_loc, 0.0_dp)
    Gp_loc = max(A_loc - Gn_loc, 0.0_dp)

    nK = size(this%kw)

    ! ---- step 2: defect self-energies, averaged over the defect line ------
    this%Sn_dc = 0.0_dp
    this%Sp_dc = 0.0_dp
    this%Sn_dv = 0.0_dp
    this%Sp_dv = 0.0_dp
    do jj = this%jlo, this%jhi
      aj = this%aline(jj - this%jlo + 1)
      do l = 0, LMAX_RUNGS
        ie_up = jj + l * this%iEshift
        ie_dn = jj - l * this%iEshift
        if (ie_up > nE .and. ie_dn < 1) exit
        if (ie_up <= nE) then
          this%Sn_dc = this%Sn_dc + aj * this%Mn(l) * Gn_loc(ie_up)
          this%Sp_dc = this%Sp_dc + aj * this%Mn(l) * this%bolt(l) * Gp_loc(ie_up)
        end if
        if (ie_dn >= 1) then
          this%Sn_dv = this%Sn_dv + aj * this%Mp(l) * this%bolt(l) * Gn_loc(ie_dn)
          this%Sp_dv = this%Sp_dv + aj * this%Mp(l) * Gp_loc(ie_dn)
        end if
      end do
    end do

    ! ---- step 3: occupation, Eq. 59 ----------------------------------------
    gam = this%Sn_dc + this%Sp_dc + this%Sn_dv + this%Sp_dv
    if (gam > 0.0_dp) then
      this%fd = (this%Sn_dc + this%Sn_dv) / gam
      this%omfd = (this%Sp_dc + this%Sp_dv) / gam
      ! ---- step 4: net recombination per defect, closed form --------------
      this%R_def = (this%Sn_dc * this%Sp_dv - this%Sp_dc * this%Sn_dv) / gam / HBAR_EV_S
    else
      this%fd = 0.5_dp
      this%omfd = 0.5_dp
      this%R_def = 0.0_dp
    end if
    this%U_rec = this%rho_d * this%R_def

    ! consistency diagnostic R_c = R_v (Eqs. 55, 58)
    Rc = (this%Sn_dc * this%omfd - this%Sp_dc * this%fd) / HBAR_EV_S
    Rv = (this%Sp_dv * this%fd - this%Sn_dv * this%omfd) / HBAR_EV_S
    gross = max(abs(this%Sn_dc * this%omfd), abs(this%Sp_dc * this%fd), &
        & abs(this%Sp_dv * this%fd), abs(this%Sn_dv * this%omfd)) / HBAR_EV_S
    if (gross > 0.0_dp) then
      bal = abs(Rc - Rv) / gross
    else
      bal = 0.0_dp
    end if

    ! ---- step 5: band self-energies from the defect Green's function ------
    xd = this%rho_d * this%cell_vol      ! defects per cell
    this%Sn_b = 0.0_dp
    this%Sp_b = 0.0_dp
    do jj = this%jlo, this%jhi
      aj = this%aline(jj - this%jlo + 1) / this%dE      ! A_d(E_j) [1/eV]
      do l = 0, LMAX_RUNGS
        ie_up = jj + l * this%iEshift      ! E = E_j + l hw : conduction-type
        ie_dn = jj - l * this%iEshift      ! E = E_j - l hw : valence-type
        if (ie_up > nE .and. ie_dn < 1) exit
        if (ie_up <= nE) then
          this%Sn_b(ie_up) = this%Sn_b(ie_up) + xd * this%Mn(l) * this%bolt(l) * this%fd * aj
          this%Sp_b(ie_up) = this%Sp_b(ie_up) + xd * this%Mn(l) * this%omfd * aj
        end if
        if (ie_dn >= 1) then
          this%Sn_b(ie_dn) = this%Sn_b(ie_dn) + xd * this%Mp(l) * this%fd * aj
          this%Sp_b(ie_dn) = this%Sp_b(ie_dn) + xd * this%Mp(l) * this%bolt(l) * this%omfd * aj
        end if
      end do
    end do

    ! ---- store Sigma^n (real positive, libNEGF convention) ----------------
    nbl = this%struct%num_PLs
    pl = this%pl_index
    orb = this%orb_index
    nrow = this%struct%mat_PL_end(pl) - this%struct%mat_PL_start(pl) + 1

    select type(pS => this%sigma_n)
    type is (TMatrixCacheMem)
      if (.not. pS%isInitialized) then
        call pS%init(nE, nK, nbl, 3, 1)
      end if
    end select

    label%spin = 1
    label%row_block = pl
    label%col_block = pl
    do iK = 1, nK
      do iE = 1, nE
        if (this%Sn_b(iE) <= 0.0_dp) cycle
        label%kpoint = iK
        label%energy_point = iE
        call create(sig, nrow, nrow)
        sig%val = (0.0_dp, 0.0_dp)
        ! same value on every orbital of the defect site
        do jj = orb, orb + this%n_orb - 1
          sig%val(jj, jj) = cmplx(this%Sn_b(iE), 0.0_dp, dp)
        end do
        call this%sigma_n%add(sig, label)
        call mat_destroy(sig)
      end do
    end do

    if (id0) then
      write(*,'(A,I3,A,F9.6,A,ES11.3,A,ES11.3,A,ES9.2)') ' SRH defect (scba_iter=', this%scba_iter, &
          & '): fd=', this%fd, '  R_def[1/s]=', this%R_def, '  U=rho_d*R[cm-3 s-1]=', this%U_rec, &
          & '  |Rc-Rv|/gross=', bal
    end if

    deallocate(Gn_loc, Gp_loc, A_loc)

  end subroutine compute_Sigma_n

  !--------------------------------------------------------------------------
  !> Sigma^r = -(i/2) (Sigma^n + Sigma^p): anti-Hermitian part fixed by the
  !! same lifetime broadening as the lesser/greater parts (current
  !! conservation); the real part (level shift) is neglected.
  !! Must be called after compute_Sigma_n of the same SCBA iteration.
  subroutine compute_Sigma_r(this, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    integer, intent(in), optional :: en_index, k_index, spin

    integer :: nE, nK, nbl, pl, orb, nrow, iE, iK
    type(TMatLabel) :: label
    type(z_dns) :: sig
    real(dp) :: gam

    if (.not. site_ok(this)) return
    if (.not. this%grid_ready) return
    if (.not. associated(this%sigma_r)) return
    if (.not. allocated(this%Sn_b)) return

    nE = this%nE
    nK = size(this%kw)
    nbl = this%struct%num_PLs
    pl = this%pl_index
    orb = this%orb_index
    nrow = this%struct%mat_PL_end(pl) - this%struct%mat_PL_start(pl) + 1

    select type(pS => this%sigma_r)
    type is (TMatrixCacheMem)
      if (.not. pS%isInitialized) then
        call pS%init(nE, nK, nbl, 3, 1)
      end if
    end select

    label%spin = 1
    label%row_block = pl
    label%col_block = pl
    do iK = 1, nK
      do iE = 1, nE
        gam = this%Sn_b(iE) + this%Sp_b(iE)
        if (gam <= 0.0_dp) cycle
        label%kpoint = iK
        label%energy_point = iE
        call create(sig, nrow, nrow)
        sig%val = (0.0_dp, 0.0_dp)
        block
          integer :: io
          do io = orb, orb + this%n_orb - 1
            sig%val(io, io) = cmplx(0.0_dp, -0.5_dp * gam, dp)
          end do
        end block
        call this%sigma_r%add(sig, label)
        call mat_destroy(sig)
      end do
    end do

  end subroutine compute_Sigma_r

  !--------------------------------------------------------------------------
  subroutine destroy_Sigma_r(this)
    class(TSrhMultiphonon) :: this
    if (associated(this%sigma_r)) then
      call this%sigma_r%destroy()
    end if
  end subroutine destroy_Sigma_r

  subroutine destroy_Sigma_n(this)
    class(TSrhMultiphonon) :: this
    if (associated(this%sigma_n)) then
      call this%sigma_n%destroy()
    end if
  end subroutine destroy_Sigma_n

  subroutine destroy(this)
    class(TSrhMultiphonon) :: this
    call this%destroy_Sigma_r()
    call this%destroy_Sigma_n()
    if (associated(this%sigma_r)) then
      deallocate(this%sigma_r)
      nullify(this%sigma_r)
    end if
    if (associated(this%sigma_n)) then
      deallocate(this%sigma_n)
      nullify(this%sigma_n)
    end if
    if (allocated(this%Mn)) deallocate(this%Mn)
    if (allocated(this%Mp)) deallocate(this%Mp)
    if (allocated(this%bolt)) deallocate(this%bolt)
    if (allocated(this%aline)) deallocate(this%aline)
    if (allocated(this%kw)) deallocate(this%kw)
    if (allocated(this%Sn_b)) deallocate(this%Sn_b)
    if (allocated(this%Sp_b)) deallocate(this%Sp_b)
  end subroutine destroy

end module srh_multiphonon
