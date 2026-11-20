!===============================================================================
! test hetero_interface
!
! PURPOSE
!   Reproduce Figs. 4, 5 and 6 of [Aeberhard2019Nonequilibrium]: spectral and
!   energy-integrated current flow at a defective AlGaAs-GaAs heterointerface
!   (about 20 nm thick, forward bias 1.1 V, same band profile in all cases):
!
!     scenario A (Fig. 4): defect plane at the interface (rho_d = 1e16 cm^-3),
!                          dark. Recombination current, expected ~ +1.8 mA/cm^2.
!     scenario B (Fig. 5): no defect, monochromatic illumination
!                          (E_gamma = 1.45 eV, I = 10 kW/m^2).
!                          Photocurrent, expected ~ -3.37 mA/cm^2.
!     scenario C (Fig. 6): defect + illumination, expected ~ -1.58 mA/cm^2
!                          (~ B + A).
!
!   For each scenario: spectral current j(E,z) [mA/(cm^2 eV)] in the conduction
!   (E in [1.2, 2.2] eV) and valence (E in [-0.2, 0.4] eV) band, and the
!   energy-integrated electron, hole and total current j_e(z), j_h(z), j_tot(z)
!   [mA/cm^2]. The total current must be constant in z (current conservation).
!   Sign convention: electric current, positive along +z.
!
! DEVICE AND BAND PROFILE (fixed input, read by eye from the figures)
!   40 cells of 0.5 nm (z = 0 ... 19.5 nm). Sites 1-10: AlGaAs, sites 11-40:
!   GaAs, interface between z = 4.5 and 5 nm. Piecewise linear band edges:
!     AlGaAs: Ec 2.08 -> 2.01, Ev 0.10 -> 0.04 eV (z = 0 -> 4.5 nm)
!     GaAs  : Ec 1.65 -> 1.45, Ev 0.23 -> 0.03 eV (z = 5 -> 19.5 nm)
!   Same effective masses on both sides (assumed): m_c = 0.067, m_v = 0.22.
!   Conduction and valence band are decoupled effective-mass chains. The
!   conduction band is contacted only at the right end (n contact, quasi-Fermi
!   level mu_c), the valence band only at the left end (p contact, mu_v), by
!   semi-infinite leads of the same chain (finite bandwidth). Inter-band
!   current can only flow via the defect or via photon absorption/emission.
!   mu_c - mu_v = 1.1 eV (bias), mu_v = Ev(0) (assumed, see below).
!
! NEGF MODEL
!   Per transverse energy and per energy E, the band chains are solved with
!   left/right-connected recursions with diagonal (local) self-energies:
!     - contact lead self-energy at the end site (complex, surface GF),
!     - intraband LO-phonon self-energies (Einstein mode, as in the p-i-n test),
!     - the defect self-energy (libnegf) at the defect cell,
!     - the photon self-energy (see below).
!   The defect Green's function is renormalised via the Dyson-Keldysh equation
!   (libnegf). SCBA iterations are carried out until convergence.
!
!   Photon self-energy (local, k-summed, absorption/stimulated emission of
!   photons of energy E_gamma in a classical field with constant coupling
!   g_ph [eV^2 cm^3]):
!     Sigma^<_c(E) = g G^<_v(E-Eg), Sigma^>_c(E) = g G^>_v(E-Eg)
!     Sigma^<_v(E) = g G^<_c(E+Eg), Sigma^>_v(E) = g G^>_c(E+Eg)
!   The paper does not give the coupling constants, so g_ph is calibrated
!   (scenario B) such that the photocurrent equals 3.37 mA/cm^2 (equivalent to
!   an absorption of 0.5 % of the incident photon flux). The same g_ph is used
!   in scenario C.
!
! SPECTRAL CURRENT
!   Due to the local scattering self-energies, the spectral current between
!   neighbouring sites follows from the continuity equation:
!     J_{n->n+1}(E) = J_{n-1->n}(E) + R_n(E),   J_{0->1} = 0,
!     R_n(E) = (1/(2 pi hbar)) sum_k [Sigma^<_n G^>_nn - Sigma^>_n G^<_nn]
!   (Sigma^<, Sigma^> total, including contacts), which is equivalent to
!   2 Re[t G^<_{n,n+1}] for the exact solution. j = -q J (electrons). The
!   current of the valence band chain is the hole current.
!
! OUTPUT
!   hetero_interface_jint.dat : z[nm]  j_e j_h j_tot (A)  j_e j_h j_tot (B)
!                               j_e j_h j_tot (C), mA/cm^2, at bond centres
!   hetero_interface_spec_{A,B,C}.dat : band(1=c,2=v)  z[nm]  E[eV]  j(E,z)
!   hetero_interface_band_edges.dat : z[nm]  Ec[eV]  Ev[eV] (cell positions)
!
! ASSUMED (not given in the reference)
!   AlGaAs masses, defect cell (first GaAs cell, i = 11), defect level
!   eps_d = Ev(z_def) + 0.75 eV, mu_v = Ev(0), photon energy rounded to the
!   energy grid (de = 3 meV -> E_gamma = 1.449 eV), T = 300 K, S = 3.5,
!   LO-phonon strength (Gamma_ph = 4 meV), eta_d = 1 meV, window +-50 meV,
!   tun_win = 50 meV, k-independent photon coupling.
!
! NOTE
!   Parameters from Fig. 3 (m_v = 0.22, eps_d = 0.75, sigma = 1e-14 cm^2) are
!   used rather than those of Table I.
!===============================================================================
program hetero_interface

  use ln_defect_srh, only : solve_defect_dyson_keldysh, &
      & compute_extended_defect_selfenergy
  use ln_defect_types, only : TDefectSpec, TDefectEnsemble, DEFECT_MODE_DELTA, &
      & DEFECT_MODE_DYSON, init_defect_spec, destruct_defect_spec, &
      & init_defect_ensemble, destruct_defect_ensemble
  use ln_defect_kernel, only : build_mpr_kernel, huang_rhys_weight
  use ln_precision, only : dp

  implicit none

  !-----------------------------------------------------------------------------
  ! Physical constants
  !-----------------------------------------------------------------------------
  real(dp), parameter :: pi = 3.14159265358979323846_dp
  real(dp), parameter :: q_e = 1.602176634e-19_dp        ! C
  real(dp), parameter :: hbar_js = 1.054571817e-34_dp    ! J s
  real(dp), parameter :: hbar_evs = 6.582119569e-16_dp   ! eV s
  real(dp), parameter :: m_e = 9.1093837015e-31_dp       ! kg
  real(dp), parameter :: kb_ev = 8.617333262e-5_dp       ! eV/K

  !-----------------------------------------------------------------------------
  ! Material and defect parameters
  !-----------------------------------------------------------------------------
  real(dp), parameter :: temp = 300.0_dp                 ! K (assumed)
  real(dp), parameter :: m_c = 0.067_dp                  ! electron mass / m0
  real(dp), parameter :: m_v = 0.22_dp                   ! hole mass / m0
  real(dp), parameter :: ec0 = 1.42_dp                   ! flat-band GaAs gap edge (calibration), eV
  real(dp), parameter :: ev0 = 0.0_dp                    ! flat-band GaAs valence edge (calibration), eV
  real(dp), parameter :: eps_d0 = 0.75_dp                ! defect level above local Ev, eV
  real(dp), parameter :: sigma_cm2 = 1.0e-14_dp          ! capture cross section, cm^2
  real(dp), parameter :: hw = 0.036_dp                   ! phonon energy, eV
  real(dp), parameter :: s_hr = 3.5_dp                   ! Huang-Rhys factor (assumed)
  real(dp), parameter :: rho_d = 1.0e16_dp               ! defect density, cm^-3

  !-----------------------------------------------------------------------------
  ! Device
  !-----------------------------------------------------------------------------
  integer,  parameter :: n_z = 40                        ! cells
  real(dp), parameter :: dz_nm = 0.5_dp                  ! cell thickness, nm
  integer,  parameter :: n_alg = 10                      ! AlGaAs cells
  integer,  parameter :: i_def = 11                      ! defect cell (first GaAs cell)
  real(dp), parameter :: vbias = 1.1_dp                  ! bias, V

  !-----------------------------------------------------------------------------
  ! Numerical parameters
  !-----------------------------------------------------------------------------
  integer,  parameter :: n_modes = 40                    ! multiphonon ladder length
  integer,  parameter :: n_en = 800                      ! energy grid points
  real(dp), parameter :: de = 0.003_dp                   ! energy step, eV (hw/de = 12)
  real(dp), parameter :: e_min = -0.2_dp                 ! lowest grid energy, eV
  integer,  parameter :: k_hw = 12                       ! hw in grid points
  integer,  parameter :: k_ph = 483                      ! photon energy in grid points (1.449 eV)
  integer,  parameter :: n_kt = 100                      ! transverse energy nodes
  real(dp), parameter :: gam = 1.0e-3_dp                 ! broadening for the calibration DOS, eV
  real(dp), parameter :: gam_lo = 4.0e-3_dp              ! LO broadening target, eV (assumed)
  integer,  parameter :: max_it = 3000                   ! SCBA iteration limit
  integer,  parameter :: min_it = 5                      ! minimum SCBA iterations
  ! The LO-phonon cascade converges slowly (ratio ~0.99 per iteration): the
  ! error of the current is ~ change/(1-0.99), and the lag between the
  ! self-energies and the Green's functions breaks local current conservation
  ! by the same amount. Hence tight tolerances for the final solutions.
  ! The energy-integrated currents are evaluated with self-energies consistent
  ! with the final Green's functions (arrays rci, rvi), so that they conserve
  ! the total current independently of the SCBA lag.
  ! The residual of the band-resolved particle conservation (sum over all cells
  ! and energies of the consistent rates rci, rvi, which vanishes exactly for a
  ! converged SCBA) is the true measure of the SCBA error and is used as
  ! convergence criterion. The self-energies are updated with Anderson mixing.
  real(dp), parameter :: tol_it = 1.0e-6_dp              ! relative tolerance on the change of the current
  real(dp), parameter :: tol_cons = 1.0e-3_dp            ! tolerance on the particle-conservation residual
  real(dp), parameter :: tol_loose = 100.0_dp            ! factor relaxing both tolerances (first calibration pass)
  integer,  parameter :: n_and = 6                       ! depth of the Anderson mixing
  real(dp), parameter :: beta_and = 0.7_dp               ! Anderson mixing parameter
  integer,  parameter :: nx = 8 * n_z * n_en + 2 * n_en  ! size of the self-energy vector
  real(dp), parameter :: eta_reg = 1.0e-10_dp            ! numerical regulator of the retarded GF, eV
  integer,  parameter :: n_omega = 41                    ! defect window points (assumed)
  real(dp), parameter :: omega_hw = 0.05_dp              ! defect window half-width, eV (assumed)
  real(dp), parameter :: eta_d = 1.0e-3_dp               ! defect broadening, eV (assumed)
  real(dp), parameter :: tun_win = 0.05_dp               ! depth of tail states coupled to the defect, eV (assumed)
  real(dp), parameter :: dmu_c = 0.0238_dp               ! Ec(n contact) - mu_c, eV (assumed; sets electron density, hence the dark current)

  !-----------------------------------------------------------------------------
  ! Photon calibration
  !-----------------------------------------------------------------------------
  real(dp), parameter :: g0_ph = 5.7e-32_dp              ! initial photon coupling, eV^2 cm^3 (from a previous calibration)
  real(dp), parameter :: j_ph_target = 3.37_dp           ! target photocurrent, mA/cm^2
  integer,  parameter :: n_cal = 8                       ! maximum number of calibration passes
  real(dp), parameter :: tol_cal = 5.0e-3_dp             ! relative accuracy of the calibrated photocurrent

  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  real(dp) :: kbt, pref_c, pref_v, neff_c, neff_v, vth_n, vth_p
  real(dp) :: cbar_n, cbar_p, avg_n, avg_p, e_l, nd_l
  real(dp) :: m0_dc, m0_dv, nb, d2_c, d2_v
  real(dp) :: dz_cm, t_c, t_v, pk_c, pk_v
  real(dp) :: ec_z(n_z), ev_z(n_z), eps_dl, ec_lo, ev_hi
  real(dp) :: mu_c, mu_v, gph, jph, r_full, cf
  real(dp) :: je_t(n_z - 1), jh_t(n_z - 1)
  real(dp) :: je_all(3, n_z - 1), jh_all(3, n_z - 1)
  real(dp) :: dev, jmean, jscale
  real(dp), allocatable :: e_grid(:), w_quad(:)
  real(dp), allocatable :: tlc(:, :), tgc(:, :), tlv(:, :), tgv(:, :)
  real(dp), allocatable :: pl_c(:, :), pg_c(:, :), pl_v(:, :), pg_v(:, :)
  real(dp), allocatable :: pp_lc(:, :), pp_gc(:, :), pp_lv(:, :), pp_gv(:, :)
  real(dp), allocatable :: rc(:, :), rv(:, :), rci(:, :), rvi(:, :)
  complex(dp), allocatable :: s_l(:, :), s_g(:, :), s_r(:, :)
  integer :: n_it, i0
  type(TDefectSpec) :: spec
  type(TDefectEnsemble) :: ens
  integer :: i, l, k, ical, isc, b, fu
  integer, parameter :: ord(3) = [2, 1, 3]
  character(len=1), parameter :: labs(3) = ['A', 'B', 'C']
  logical :: sdef, sph

  !=============================================================================
  ! STEP 1: Constants and band quantities
  !=============================================================================
  kbt = kb_ev * temp
  pref_c = dos_pref(m_c)
  pref_v = dos_pref(m_v)
  neff_c = pref_c * 0.5_dp * sqrt(pi) * kbt**1.5_dp      ! effective DOS, cm^-3
  neff_v = pref_v * 0.5_dp * sqrt(pi) * kbt**1.5_dp
  vth_n = 100.0_dp * sqrt(3.0_dp * kbt * q_e / (m_c * m_e))   ! cm/s
  vth_p = 100.0_dp * sqrt(3.0_dp * kbt * q_e / (m_v * m_e))
  cbar_n = vth_n * sigma_cm2                             ! cm^3/s
  cbar_p = vth_p * sigma_cm2

  dz_cm = dz_nm * 1.0e-7_dp
  t_c = hbar_js**2 / (2.0_dp * m_c * m_e * (dz_nm * 1.0e-9_dp)**2) / q_e
  t_v = hbar_js**2 / (2.0_dp * m_v * m_e * (dz_nm * 1.0e-9_dp)**2) / q_e
  pk_c = m_c * m_e / (pi * hbar_js**2) * q_e * 1.0e-4_dp
  pk_v = m_v * m_e / (pi * hbar_js**2) * q_e * 1.0e-4_dp

  !=============================================================================
  ! STEP 2: Calibrate the multiphonon prefactor (flat-band GaAs bulk)
  !=============================================================================
  call init_defect_spec(spec, 1, [1], [1.0_dp], n_modes)
  call configure_spec(spec, 1.0_dp, 1.0_dp, 1.0_dp)
  avg_n = 0.0_dp
  avg_p = 0.0_dp
  do l = 1, n_modes
    e_l = eps_d0 + spec%e_shift(l)
    if (e_l >= ec0) then
      nd_l = chain_dos(1, e_l - ec0)
      avg_n = avg_n + 2.0_dp * pi / hbar_evs * spec%m_dc_capt(l) * nd_l &
          & * exp(-(e_l - ec0) / kbt)
    end if
    e_l = eps_d0 - spec%e_shift(l)
    if (e_l <= ev0) then
      nd_l = chain_dos(2, ev0 - e_l)
      avg_p = avg_p + 2.0_dp * pi / hbar_evs * spec%m_dv_em(l) * nd_l &
          & * exp(-(ev0 - e_l) / kbt)
    end if
  end do
  avg_n = avg_n / neff_c
  avg_p = avg_p / neff_v
  m0_dc = cbar_n / avg_n
  m0_dv = cbar_p / avg_p
  call destruct_defect_spec(spec)

  ! LO-phonon coupling of the bands (Einstein mode)
  nb = 1.0_dp / (exp(hw / kbt) - 1.0_dp)
  d2_c = gam_lo / ((nb + 1.0_dp) * 2.0_dp * pi * pref_c * sqrt(hw))
  d2_v = gam_lo / ((nb + 1.0_dp) * 2.0_dp * pi * pref_v * sqrt(hw))

  write(*, '(a)') '--- parameters ------------------------------------------'
  write(*, '(a,es12.4,a,es12.4)') ' v_th*sigma (n,p) [cm^3/s] = ', cbar_n, &
      & '   ', cbar_p
  write(*, '(a,es12.4,a,es12.4)') ' M0_dc = ', m0_dc, '   M0_dv = ', m0_dv
  write(*, '(a,es12.4,a,es12.4)') ' t_c [eV] = ', t_c, '   t_v [eV] = ', t_v

  !=============================================================================
  ! STEP 3: Band profile, energy grid, defect ensemble, arrays
  !=============================================================================
  do i = 1, n_alg
    ec_z(i) = 2.08_dp - 0.07_dp * real(i - 1, dp) / real(n_alg - 1, dp)
    ev_z(i) = 0.10_dp - 0.06_dp * real(i - 1, dp) / real(n_alg - 1, dp)
  end do
  do i = n_alg + 1, n_z
    ec_z(i) = 1.65_dp - 0.20_dp * real(i - n_alg - 1, dp) / real(n_z - n_alg - 1, dp)
    ev_z(i) = 0.23_dp - 0.20_dp * real(i - n_alg - 1, dp) / real(n_z - n_alg - 1, dp)
  end do

  mu_c = ec_z(n_z) - dmu_c
  mu_v = mu_c - vbias

  eps_dl = ev_z(i_def) + eps_d0
  i0 = nint((eps_dl - e_min) / de) + 1
  ec_lo = minval(ec_z) - 0.1_dp
  ev_hi = maxval(ev_z) + 0.1_dp

  call init_defect_ensemble(ens, 1, 1)
  call init_defect_spec(ens%defects(1), 1, [1], [1.0_dp], n_modes)
  call configure_spec(ens%defects(1), m0_dc, m0_dv, rho_d)
  ens%defects(1)%mode = DEFECT_MODE_DYSON
  ens%defects(1)%n_omega = n_omega
  ens%defects(1)%omega_halfwidth = omega_hw
  ens%defects(1)%eta_d = eta_d
  ens%defects(1)%eps_d = eps_dl
  ens%defects(1)%ec_cut = ec_z(i_def) - tun_win
  ens%defects(1)%ev_cut = ev_z(i_def) + tun_win

  allocate(e_grid(n_en), w_quad(n_en))
  allocate(tlc(n_z, n_en), tgc(n_z, n_en), tlv(n_z, n_en), tgv(n_z, n_en))
  allocate(pl_c(n_z, n_en), pg_c(n_z, n_en), pl_v(n_z, n_en), pg_v(n_z, n_en))
  allocate(pp_lc(n_z, n_en), pp_gc(n_z, n_en), pp_lv(n_z, n_en), pp_gv(n_z, n_en))
  allocate(rc(n_z, n_en), rv(n_z, n_en), rci(n_z, n_en), rvi(n_z, n_en))
  allocate(s_l(1, n_en), s_g(1, n_en), s_r(1, n_en))
  w_quad(:) = de
  do i = 1, n_en
    e_grid(i) = eps_dl + de * real(i - i0, dp)
  end do

  cf = -q_e / (2.0_dp * pi * hbar_evs) * 1.0e3_dp * de   ! mA/cm^2 per cm^-2

  !=============================================================================
  ! STEP 4: Scenarios B (calibration of g_ph), A, C
  !=============================================================================
  gph = g0_ph
  do k = 1, 3
    isc = ord(k)
    sdef = (isc /= 2)
    sph = (isc /= 1)

    if (isc == 2) then
      ! photocurrent is linear in g_ph: iterate until the tightly converged
      ! pass reproduces the target (the first pass is loose, only a guess)
      do ical = 1, n_cal
        call negf_scba(sdef, sph, gph, (ical > 1), r_full, n_it)
        call currents(je_t, jh_t)
        jph = je_t(n_z - 1) + jh_t(n_z - 1)
        write(*, '(a,i0,a,es12.4,a,es14.6,a,i0)') ' calibration pass ', ical, &
            & '  g_ph=', gph, '  j_ph=', jph, '  n_iter=', n_it
        if (ical > 1 .and. abs(abs(jph) - j_ph_target) <= tol_cal * j_ph_target) exit
        if (ical < n_cal) then
          if (abs(jph) > 1.0e-30_dp) then
            gph = gph * j_ph_target / abs(jph)
          else
            gph = gph * 10.0_dp
          end if
        end if
      end do
    else
      call negf_scba(sdef, sph, gph, .true., r_full, n_it)
      call currents(je_t, jh_t)
    end if

    je_all(isc, :) = je_t(:)
    jh_all(isc, :) = jh_t(:)
    write(*, '(a,a,a,es14.6,a,es14.6,a,i0)') ' scenario ', labs(isc), &
        & ': j_tot(1)=', je_t(1) + jh_t(1), '  j_tot(end)=', &
        & je_t(n_z - 1) + jh_t(n_z - 1), '  n_iter=', n_it
    call write_spectrum(labs(isc))
  end do

  call destruct_defect_ensemble(ens)

  !=============================================================================
  ! STEP 5: Output
  !=============================================================================
  open(newunit=fu, file='hetero_interface_jint.dat', status='replace', action='write')
  write(fu, '(a)') '# Energy-integrated currents, mA/cm^2, at bond centres'
  write(fu, '(a)') '# z[nm] | A: j_e j_h j_tot | B: j_e j_h j_tot | C: j_e j_h j_tot'
  do b = 1, n_z - 1
    write(fu, '(f8.3,9es16.7e3)') (real(b, dp) - 0.5_dp) * dz_nm, &
        & je_all(1, b), jh_all(1, b), je_all(1, b) + jh_all(1, b), &
        & je_all(2, b), jh_all(2, b), je_all(2, b) + jh_all(2, b), &
        & je_all(3, b), jh_all(3, b), je_all(3, b) + jh_all(3, b)
  end do
  close(fu)

  ! Export the exact band-edge profile used by the transport Hamiltonians.
  ! Positions are the longitudinal site coordinates (0, 0.5, ..., 19.5 nm).
  open(newunit=fu, file='hetero_interface_band_edges.dat', status='replace', action='write')
  write(fu, '(a)') '# z[nm]  Ec[eV]  Ev[eV]'
  do i = 1, n_z
    write(fu, '(f8.3,2f12.6)') real(i - 1, dp) * dz_nm, ec_z(i), ev_z(i)
  end do
  close(fu)

  write(*, '(a)') '--- results (mA/cm^2) -----------------------------------'
  write(*, '(a)') ' scen   j_tot(mean)   target    j_e(first)  j_h(first)'
  write(*, '(2x,a,es14.5,a,es14.5,es14.5)') 'A', mean_jtot(1), '   +1.8 ', &
      & je_all(1, 1), jh_all(1, 1)
  write(*, '(2x,a,es14.5,a,es14.5,es14.5)') 'B', mean_jtot(2), '   -3.37', &
      & je_all(2, 1), jh_all(2, 1)
  write(*, '(2x,a,es14.5,a,es14.5,es14.5)') 'C', mean_jtot(3), '   -1.58', &
      & je_all(3, 1), jh_all(3, 1)

  !=============================================================================
  ! STEP 6: Checks
  !=============================================================================
  do isc = 1, 3
    jmean = mean_jtot(isc)
    jscale = max(maxval(abs(je_all(isc, :))), maxval(abs(jh_all(isc, :))), 1.0e-30_dp)
    dev = 0.0_dp
    do b = 1, n_z - 1
      dev = max(dev, abs(je_all(isc, b) + jh_all(isc, b) - jmean))
    end do
    if (dev > 5.0e-3_dp * jscale) then
      write(*, '(a,a,a,es10.2)') 'FAIL: total current not conserved in scenario ', &
          & labs(isc), ', relative deviation ', dev / jscale
      error stop 1
    end if
  end do
  if (mean_jtot(1) <= 0.0_dp) then
    write(*, '(a)') 'FAIL: dark recombination current (A) not positive'
    error stop 1
  end if
  if (mean_jtot(2) >= 0.0_dp) then
    write(*, '(a)') 'FAIL: photocurrent (B) not negative'
    error stop 1
  end if
  if (mean_jtot(3) <= mean_jtot(2)) then
    write(*, '(a)') 'FAIL: current of (C) not above the pure photocurrent (B)'
    error stop 1
  end if

  write(*, '(a)') 'test f90hetero_interface finished, data in hetero_interface_*.dat'

contains

  !-----------------------------------------------------------------------------
  ! Mean total current of scenario s over the bonds
  !-----------------------------------------------------------------------------
  function mean_jtot(s) result(m)
    integer, intent(in) :: s
    real(dp) :: m
    m = sum(je_all(s, :) + jh_all(s, :)) / real(n_z - 1, dp)
  end function mean_jtot

  !-----------------------------------------------------------------------------
  ! Fermi-Dirac occupation of x = (E - mu)/kT, overflow-safe
  !-----------------------------------------------------------------------------
  pure function fermi(x) result(f)
    real(dp), intent(in) :: x
    real(dp) :: f
    if (x > 500.0_dp) then
      f = 0.0_dp
    else if (x < -500.0_dp) then
      f = 1.0_dp
    else
      f = 1.0_dp / (1.0_dp + exp(x))
    end if
  end function fermi

  !-----------------------------------------------------------------------------
  ! 3D parabolic density of states N(E) = pref * sqrt(E), [1/(eV cm^3)]
  !-----------------------------------------------------------------------------
  pure function dos_pref(m_rel) result(pref)
    real(dp), intent(in) :: m_rel
    real(dp) :: pref
    pref = (2.0_dp * m_rel * m_e * q_e / hbar_js**2)**1.5_dp &
        & / (2.0_dp * pi**2) * 1.0e-6_dp
  end function dos_pref

  !-----------------------------------------------------------------------------
  ! DOS [1/(eV cm^3)] of the discrete flat-band chain with transverse midpoint
  ! quadrature and broadening gam (calibration only). ib = 1: conduction,
  ! ib = 2: valence.
  !-----------------------------------------------------------------------------
  pure function chain_dos(ib, x_in) result(nd)
    integer, intent(in) :: ib
    real(dp), intent(in) :: x_in
    real(dp) :: nd
    real(dp) :: t_h, pk_h, x, eperp
    complex(dp) :: zeta, gr
    integer :: j

    if (ib == 1) then
      t_h = t_c
      pk_h = pk_c
    else
      t_h = t_v
      pk_h = pk_v
    end if
    nd = 0.0_dp
    do j = 1, n_kt
      eperp = (real(j, dp) - 0.5_dp) * de
      x = x_in - eperp
      if (x < 0.0_dp) exit
      if (x > 4.0_dp * t_h) cycle
      zeta = cmplx(x - 2.0_dp * t_h, 0.5_dp * gam, dp)
      gr = 1.0_dp / sqrt(zeta * zeta - 4.0_dp * t_h**2)
      if (aimag(gr) > 0.0_dp) gr = -gr
      nd = nd - aimag(gr) / pi
    end do
    nd = nd * de * pk_h / dz_cm
  end function chain_dos

  !-----------------------------------------------------------------------------
  ! Multiphonon line-shape weights, Schenk Eq. (12) with r_F = 0
  !-----------------------------------------------------------------------------
  pure function shape_n(l) result(s)
    integer, intent(in) :: l
    real(dp) :: s
    s = (real(l, dp) - s_hr)**2 / s_hr * huang_rhys_weight(l, s_hr, hw, kbt)
  end function shape_n

  pure function shape_p(l) result(s)
    integer, intent(in) :: l
    real(dp) :: s
    s = (real(l, dp) + s_hr)**2 / s_hr * huang_rhys_weight(l, s_hr, hw, kbt)
  end function shape_p

  !-----------------------------------------------------------------------------
  ! Set the defect parameters and build the libnegf MPR kernel
  !-----------------------------------------------------------------------------
  subroutine configure_spec(sp, m0c, m0v, dens)
    type(TDefectSpec), intent(inout) :: sp
    real(dp), intent(in) :: m0c, m0v, dens
    sp%eps_d = eps_d0
    sp%kbt = kbt
    sp%density_weight = dens
    sp%hbar_omega0 = hw
    sp%huang_rhys_n = s_hr
    sp%huang_rhys_p = s_hr
    sp%rf_rph_ratio_n = 0.0_dp
    sp%rf_rph_ratio_p = 0.0_dp
    sp%m0_dc = m0c
    sp%m0_dv = m0v
    sp%ec_cut = ec0
    sp%ev_cut = ev0
    sp%mode = DEFECT_MODE_DELTA
    sp%pl_index = 1
    call build_mpr_kernel(sp)
  end subroutine configure_spec

  !-----------------------------------------------------------------------------
  ! Retarded self-energy of a semi-infinite lead (same chain, hopping t_h) at the
  ! end site; z_in = E_z - h_site (h_site includes the +-2t band offset).
  !-----------------------------------------------------------------------------
  subroutine contact_se(t_h, z_in, sig)
    real(dp), intent(in) :: t_h, z_in
    complex(dp), intent(out) :: sig

    complex(dp) :: zc, sq, s1, s2

    zc = cmplx(z_in, 1.0e-6_dp, dp)
    sq = sqrt(zc * zc - 4.0_dp * t_h**2)
    s1 = 0.5_dp * (zc - sq)
    s2 = 0.5_dp * (zc + sq)
    if (abs(abs(s1) - abs(s2)) > 1.0e-6_dp * t_h) then
      if (abs(s1) < abs(s2)) then
        sig = s1
      else
        sig = s2
      end if
    else
      if (aimag(s1) < aimag(s2)) then
        sig = s1
      else
        sig = s2
      end if
    end if
  end subroutine contact_se

  !-----------------------------------------------------------------------------
  ! Transverse-momentum summed local Green's functions of one band at energy
  ! index ie: tilde G<(z), tilde G>(z) in cm^-3 eV^-1 (G< = i tilde G<,
  ! G> = -i tilde G>), and the spectral net in-scattering rate per area
  ! rro(z) = sum_k [Sigma~^< G~^> - Sigma~^> G~^<] [cm^-2].
  !   ib = 1: conduction band chain (contact at site n_z),
  !   ib = 2: valence band chain (contact at site 1).
  !-----------------------------------------------------------------------------
  subroutine band_energy(ib, ie, mu, tlo, tgo, rro)
    integer, intent(in) :: ib, ie
    real(dp), intent(in) :: mu
    real(dp), intent(out) :: tlo(n_z), tgo(n_z), rro(n_z)

    real(dp) :: ee, fo, fh, t2, t_h, wk, ka, ez, eperp, w_site, gamc, gl_n, gg_n
    real(dp) :: hh(n_z), tsl(n_z), tsg(n_z), ag(n_z)
    real(dp) :: tslj(n_z), tsgj(n_z)
    real(dp) :: laml(n_z), rhol(n_z), lamg(n_z), rhog(n_z)
    complex(dp) :: sg(n_z), glc(n_z), grc(n_z), den, sig
    integer :: j, n, nc

    ee = e_grid(ie)
    fo = fermi((ee - mu) / kbt)
    fh = 1.0_dp - fo

    if (ib == 1) then
      hh(:) = ec_z(:) + 2.0_dp * t_c
      t_h = t_c
      ka = pk_c * de
      tsl(:) = pl_c(:, ie) + pp_lc(:, ie)
      tsg(:) = pg_c(:, ie) + pp_gc(:, ie)
      nc = n_z
    else
      hh(:) = ev_z(:) - 2.0_dp * t_v
      t_h = t_v
      ka = pk_v * de
      tsl(:) = pl_v(:, ie) + pp_lv(:, ie)
      tsg(:) = pg_v(:, ie) + pp_gv(:, ie)
      nc = 1
    end if
    t2 = t_h**2
    wk = ka / dz_cm

    ! defect self-energy at the defect cell (libnegf). libnegf includes the
    ! factor 1/(2 pi) of the omega integral (G_d normalised to 1 with dw/2pi)
    ! in sigma_ext, i.e. sigma_ext = Sigma_paper/(2 pi); the paper's Eqs. 18-19
    ! have no 1/(2 pi), hence the factor 2 pi.
    if ((ib == 1 .and. ee > eps_dl) .or. (ib == 2 .and. ee < eps_dl)) then
      tsl(i_def) = tsl(i_def) + 2.0_dp * pi * aimag(s_l(1, ie))
      tsg(i_def) = tsg(i_def) - 2.0_dp * pi * aimag(s_g(1, ie))
    end if

    tlo(:) = 0.0_dp
    tgo(:) = 0.0_dp
    rro(:) = 0.0_dp

    do j = 1, n_kt
      eperp = (real(j, dp) - 0.5_dp) * de
      if (ib == 1) then
        ez = ee - eperp
      else
        ez = ee + eperp
      end if

      do n = 1, n_z
        tslj(n) = tsl(n)
        tsgj(n) = tsg(n)
        w_site = 0.5_dp * (tslj(n) + tsgj(n)) + eta_reg
        sg(n) = cmplx(0.0_dp, -w_site, dp)
      end do

      ! contact lead at the end site
      call contact_se(t_h, ez - hh(nc), sig)
      gamc = max(-2.0_dp * aimag(sig), 0.0_dp)
      sg(nc) = sg(nc) + sig
      tslj(nc) = tslj(nc) + gamc * fo
      tsgj(nc) = tsgj(nc) + gamc * fh

      ! left- and right-connected retarded functions
      glc(1) = 1.0_dp / (ez - hh(1) - sg(1))
      do n = 2, n_z
        glc(n) = 1.0_dp / (ez - hh(n) - sg(n) - t2 * glc(n - 1))
      end do
      grc(n_z) = 1.0_dp / (ez - hh(n_z) - sg(n_z))
      do n = n_z - 1, 1, -1
        grc(n) = 1.0_dp / (ez - hh(n) - sg(n) - t2 * grc(n + 1))
      end do

      ! |G_nn|^2
      do n = 1, n_z
        den = ez - hh(n) - sg(n)
        if (n > 1) den = den - t2 * glc(n - 1)
        if (n < n_z) den = den - t2 * grc(n + 1)
        ag(n) = 1.0_dp / (real(den, dp)**2 + aimag(den)**2)
      end do

      ! sums over the left and right sites of |G_nm|^2 * tilde Sigma_m
      laml(1) = tslj(1)
      lamg(1) = tsgj(1)
      do n = 2, n_z
        laml(n) = tslj(n) + t2 * (real(glc(n - 1), dp)**2 &
            & + aimag(glc(n - 1))**2) * laml(n - 1)
        lamg(n) = tsgj(n) + t2 * (real(glc(n - 1), dp)**2 &
            & + aimag(glc(n - 1))**2) * lamg(n - 1)
      end do
      rhol(n_z) = tslj(n_z)
      rhog(n_z) = tsgj(n_z)
      do n = n_z - 1, 1, -1
        rhol(n) = tslj(n) + t2 * (real(grc(n + 1), dp)**2 &
            & + aimag(grc(n + 1))**2) * rhol(n + 1)
        rhog(n) = tsgj(n) + t2 * (real(grc(n + 1), dp)**2 &
            & + aimag(grc(n + 1))**2) * rhog(n + 1)
      end do

      do n = 1, n_z
        gl_n = ag(n) * (laml(n) + rhol(n) - tslj(n))
        gg_n = ag(n) * (lamg(n) + rhog(n) - tsgj(n))
        tlo(n) = tlo(n) + wk * gl_n
        tgo(n) = tgo(n) + wk * gg_n
        rro(n) = rro(n) + ka * (tslj(n) * gg_n - tsgj(n) * gl_n)
      end do
    end do
  end subroutine band_energy

  !-----------------------------------------------------------------------------
  ! NEGF self-consistent (SCBA) cycle of the heterostructure.
  !   sdef: defect on, sph: photon coupling on with strength gp.
  ! On return: r_conv = |J_e(n_z-1)| + |J_v(1)| (energy-integrated particle
  ! fluxes at the contacts, cm^-2), nit = passes used. rc, rv hold the spectral
  ! net in-scattering rates of the last pass.
  !-----------------------------------------------------------------------------
  subroutine negf_scba(sdef_in, sph_in, gp, tight, r_conv, nit)
    logical, intent(in) :: sdef_in, sph_in, tight
    real(dp), intent(in) :: gp
    real(dp), intent(out) :: r_conv
    integer, intent(out) :: nit

    complex(dp) :: gl(1, n_en), gg(1, n_en)
    complex(dp) :: sl_r(1, n_en), sg_r(1, n_en)
    real(dp) :: nlc(n_z), ngc(n_z), nlv(n_z), ngv(n_z)
    real(dp) :: ee, r_old, r_new, jt_first, jt_last, mism, tl_it, tl_cons
    real(dp) :: dl_c, dg_c, dl_v, dg_v
    real(dp), allocatable :: xo(:), xg(:), xn(:), fcur(:), xprev(:), fprev(:)
    real(dp), allocatable :: dxh(:, :), dfh(:, :)
    real(dp) :: amat(n_and, n_and), bvec(n_and), gam_a(n_and), cmax
    integer :: ie, it, ip, im, nn, nh, ih, jh, mm
    logical :: have_prev, conv
    if (tight) then
      tl_it = tol_it
      tl_cons = tol_cons
    else
      tl_it = tol_loose * tol_it
      tl_cons = tol_loose * tol_cons
    end if

    allocate(xo(nx), xg(nx), xn(nx), fcur(nx), xprev(nx), fprev(nx))
    allocate(dxh(nx, n_and), dfh(nx, n_and))
    nh = 0
    have_prev = .false.
    s_l(:, :) = (0.0_dp, 0.0_dp)
    s_g(:, :) = (0.0_dp, 0.0_dp)
    s_r(:, :) = (0.0_dp, 0.0_dp)
    pl_c(:, :) = 0.0_dp
    pg_c(:, :) = 0.0_dp
    pl_v(:, :) = 0.0_dp
    pg_v(:, :) = 0.0_dp
    pp_lc(:, :) = 0.0_dp
    pp_gc(:, :) = 0.0_dp
    pp_lv(:, :) = 0.0_dp
    pp_gv(:, :) = 0.0_dp
    r_old = -1.0_dp
    r_new = 0.0_dp
    nit = max_it
    conv = .false.

    do it = 1, max_it
      ! -- band Green's functions with the current self-energies ------------
      tlc(:, :) = 0.0_dp
      tgc(:, :) = 0.0_dp
      tlv(:, :) = 0.0_dp
      tgv(:, :) = 0.0_dp
      rc(:, :) = 0.0_dp
      rv(:, :) = 0.0_dp
      do ie = 1, n_en
        ee = e_grid(ie)
        if (ee >= ec_lo) call band_energy(1, ie, mu_c, tlc(:, ie), tgc(:, ie), rc(:, ie))
        if (ee <= ev_hi) call band_energy(2, ie, mu_v, tlv(:, ie), tgv(:, ie), rv(:, ie))
        gl(1, ie) = cmplx(0.0_dp, tlc(i_def, ie) + tlv(i_def, ie), dp)
        gg(1, ie) = cmplx(0.0_dp, -(tgc(i_def, ie) + tgv(i_def, ie)), dp)
      end do

      ! -- energy-integrated rates: remove the scattering part evaluated with
      !    the old self-energies (it is re-added below with the new ones) ----
      do ie = 1, n_en
        ee = e_grid(ie)
        do nn = 1, n_z
          rci(nn, ie) = rc(nn, ie) - dz_cm * ((pl_c(nn, ie) + pp_lc(nn, ie)) * tgc(nn, ie) &
              & - (pg_c(nn, ie) + pp_gc(nn, ie)) * tlc(nn, ie))
          rvi(nn, ie) = rv(nn, ie) - dz_cm * ((pl_v(nn, ie) + pp_lv(nn, ie)) * tgv(nn, ie) &
              & - (pg_v(nn, ie) + pp_gv(nn, ie)) * tlv(nn, ie))
        end do
        if (ee > eps_dl) then
          dl_c = 2.0_dp * pi * aimag(s_l(1, ie))
          dg_c = -2.0_dp * pi * aimag(s_g(1, ie))
          rci(i_def, ie) = rci(i_def, ie) - dz_cm * (dl_c * tgc(i_def, ie) - dg_c * tlc(i_def, ie))
        end if
        if (ee < eps_dl) then
          dl_v = 2.0_dp * pi * aimag(s_l(1, ie))
          dg_v = -2.0_dp * pi * aimag(s_g(1, ie))
          rvi(i_def, ie) = rvi(i_def, ie) - dz_cm * (dl_v * tgv(i_def, ie) - dg_v * tlv(i_def, ie))
        end if
      end do

      ! -- convergence quantity: energy-integrated fluxes at the contacts ----
      r_new = abs(sum(rc(1:n_z - 1, :))) * de + abs(sum(rv(1, :))) * de

      ! -- old self-energy vector (input of this pass) ------------------------
      call pack_x(xo)

      ! -- defect occupation and self-energy (libnegf) -----------------------
      s_l(:, :) = (0.0_dp, 0.0_dp)
      s_g(:, :) = (0.0_dp, 0.0_dp)
      s_r(:, :) = (0.0_dp, 0.0_dp)
      sl_r(:, :) = (0.0_dp, 0.0_dp)
      sg_r(:, :) = (0.0_dp, 0.0_dp)
      if (sdef_in) then
        call solve_defect_dyson_keldysh(ens%defects(1), e_grid, w_quad, gl, gg)
        call compute_extended_defect_selfenergy(ens, e_grid, w_quad, s_l, s_g, s_r)
        ! unmixed self-energy, consistent with the defect balance of this pass
        sl_r(:, :) = s_l(:, :)
        sg_r(:, :) = s_g(:, :)
      end if

      ! -- LO-phonon and photon self-energies from the local G (mixed) -------
      do ie = 1, n_en
        ip = ie + k_hw
        im = ie - k_hw
        nlc(:) = 0.0_dp
        ngc(:) = 0.0_dp
        nlv(:) = 0.0_dp
        ngv(:) = 0.0_dp
        if (ip <= n_en) then
          nlc(:) = nlc(:) + (nb + 1.0_dp) * tlc(:, ip)
          ngc(:) = ngc(:) + nb * tgc(:, ip)
          nlv(:) = nlv(:) + (nb + 1.0_dp) * tlv(:, ip)
          ngv(:) = ngv(:) + nb * tgv(:, ip)
        end if
        if (im >= 1) then
          nlc(:) = nlc(:) + nb * tlc(:, im)
          ngc(:) = ngc(:) + (nb + 1.0_dp) * tgc(:, im)
          nlv(:) = nlv(:) + nb * tlv(:, im)
          ngv(:) = ngv(:) + (nb + 1.0_dp) * tgv(:, im)
        end if
        pl_c(:, ie) = d2_c * nlc(:)
        pg_c(:, ie) = d2_c * ngc(:)
        pl_v(:, ie) = d2_v * nlv(:)
        pg_v(:, ie) = d2_v * ngv(:)

        if (sph_in) then
          nlc(:) = 0.0_dp
          ngc(:) = 0.0_dp
          nlv(:) = 0.0_dp
          ngv(:) = 0.0_dp
          if (ie - k_ph >= 1) then
            nlc(:) = gp * tlv(:, ie - k_ph)
            ngc(:) = gp * tgv(:, ie - k_ph)
          end if
          if (ie + k_ph <= n_en) then
            nlv(:) = gp * tlc(:, ie + k_ph)
            ngv(:) = gp * tgc(:, ie + k_ph)
          end if
          ! no absorption in the AlGaAs cells (E_gamma < gap there): the
          ! photocurrent must vanish exactly for z < interface (Figs. 5, 6)
          nlc(1:n_alg) = 0.0_dp
          ngc(1:n_alg) = 0.0_dp
          nlv(1:n_alg) = 0.0_dp
          ngv(1:n_alg) = 0.0_dp
          pp_lc(:, ie) = nlc(:)
          pp_gc(:, ie) = ngc(:)
          pp_lv(:, ie) = nlv(:)
          pp_gv(:, ie) = ngv(:)
        end if
      end do

      ! -- energy-integrated rates: re-add the scattering part with the new
      !    self-energies (consistent with the Green's functions of this pass).
      !    The defect part uses the unmixed self-energy (sl_r, sg_r).
      do ie = 1, n_en
        ee = e_grid(ie)
        do nn = 1, n_z
          rci(nn, ie) = rci(nn, ie) + dz_cm * ((pl_c(nn, ie) + pp_lc(nn, ie)) * tgc(nn, ie) &
              & - (pg_c(nn, ie) + pp_gc(nn, ie)) * tlc(nn, ie))
          rvi(nn, ie) = rvi(nn, ie) + dz_cm * ((pl_v(nn, ie) + pp_lv(nn, ie)) * tgv(nn, ie) &
              & - (pg_v(nn, ie) + pp_gv(nn, ie)) * tlv(nn, ie))
        end do
        if (ee > eps_dl) then
          dl_c = 2.0_dp * pi * aimag(sl_r(1, ie))
          dg_c = -2.0_dp * pi * aimag(sg_r(1, ie))
          rci(i_def, ie) = rci(i_def, ie) + dz_cm * (dl_c * tgc(i_def, ie) - dg_c * tlc(i_def, ie))
        end if
        if (ee < eps_dl) then
          dl_v = 2.0_dp * pi * aimag(sl_r(1, ie))
          dg_v = -2.0_dp * pi * aimag(sg_r(1, ie))
          rvi(i_def, ie) = rvi(i_def, ie) + dz_cm * (dl_v * tgv(i_def, ie) - dg_v * tlv(i_def, ie))
        end if
      end do

      ! -- convergence: small change AND current conservation --------------
      ! (the change per iteration alone is misleading for slowly converging
      ! LO-phonon cascades; the mismatch between the total current at the
      ! first and last bond measures the true residual error)
      ! sum over all cells and energies of rci (rvi) = sum of (new - old)
      ! scattering rates of that band: exactly zero for a converged SCBA
      ! (it equals the spurious hole/electron current leaking past the closed
      ! end of a band chain, e.g. j_h at z > z_defect in scenario A)
      mism = (abs(sum(rci)) + abs(sum(rvi))) * de / (r_new + 1.0e-30_dp)
      if (mod(it, 20) == 0) then
        write(*, '(a,i0,a,es14.6,a,es10.2,a,es10.2)') '   scba it=', it, &
            & '  flux=', r_new, '  change=', &
            & abs(r_new - r_old) / (abs(r_new) + 1.0e-30_dp), '  mismatch=', mism
      end if
      if (it >= min_it .and. r_old >= 0.0_dp) then
        if (abs(r_new - r_old) <= tl_it * abs(r_new) + 1.0e-30_dp &
            & .and. mism <= tl_cons) then
          nit = it
          exit
        end if
      end if
      ! -- Anderson mixing of the self-energy vector -------------------------
      !    x_new = x_old + beta f - (dX + beta dF) gamma,  f = g(x_old) - x_old
      if (it < max_it) then
        call pack_x(xg)
        fcur(:) = xg(:) - xo(:)
        if (have_prev) then
          ih = mod(nh, n_and) + 1
          dxh(:, ih) = xo(:) - xprev(:)
          dfh(:, ih) = fcur(:) - fprev(:)
          nh = nh + 1
        end if
        xprev(:) = xo(:)
        fprev(:) = fcur(:)
        have_prev = .true.
        mm = min(nh, n_and)
        xn(:) = xo(:) + beta_and * fcur(:)
        if (mm > 0) then
          do ih = 1, mm
            do jh = 1, ih
              amat(ih, jh) = dot_product(dfh(:, ih), dfh(:, jh))
              amat(jh, ih) = amat(ih, jh)
            end do
            bvec(ih) = dot_product(dfh(:, ih), fcur(:))
          end do
          cmax = 0.0_dp
          do ih = 1, mm
            cmax = max(cmax, amat(ih, ih))
          end do
          do ih = 1, mm
            amat(ih, ih) = amat(ih, ih) + 1.0e-10_dp * cmax + 1.0e-300_dp
          end do
          call solve_small(mm, amat, bvec, gam_a)
          do ih = 1, mm
            xn(:) = xn(:) - gam_a(ih) * (dxh(:, ih) + beta_and * dfh(:, ih))
          end do
        end if
        call unpack_x(xn)
      end if
      r_old = r_new
    end do

    deallocate(xo, xg, xn, fcur, xprev, fprev, dxh, dfh)

    if (nit >= max_it) then
      write(*, '(a,es10.2)') 'WARNING: SCBA not converged, last relative change ', &
          & abs(r_new - r_old) / (abs(r_new) + 1.0e-30_dp)
    end if

    r_conv = r_new
  end subroutine negf_scba

  !-----------------------------------------------------------------------------
  ! Pack / unpack the self-energies (phonon, photon, defect) into one vector.
  ! The defect part stores aimag(s_l), aimag(s_g) (only these are used).
  !-----------------------------------------------------------------------------
  subroutine pack_x(v)
    real(dp), intent(out) :: v(nx)
    integer :: o, nze

    nze = n_z * n_en
    o = 0
    v(o + 1:o + nze) = reshape(pl_c, [nze]); o = o + nze
    v(o + 1:o + nze) = reshape(pg_c, [nze]); o = o + nze
    v(o + 1:o + nze) = reshape(pl_v, [nze]); o = o + nze
    v(o + 1:o + nze) = reshape(pg_v, [nze]); o = o + nze
    v(o + 1:o + nze) = reshape(pp_lc, [nze]); o = o + nze
    v(o + 1:o + nze) = reshape(pp_gc, [nze]); o = o + nze
    v(o + 1:o + nze) = reshape(pp_lv, [nze]); o = o + nze
    v(o + 1:o + nze) = reshape(pp_gv, [nze]); o = o + nze
    v(o + 1:o + n_en) = aimag(s_l(1, :)); o = o + n_en
    v(o + 1:o + n_en) = aimag(s_g(1, :))
  end subroutine pack_x

  subroutine unpack_x(v)
    real(dp), intent(in) :: v(nx)
    integer :: o, nze

    nze = n_z * n_en
    o = 0
    pl_c = max(reshape(v(o + 1:o + nze), [n_z, n_en]), 0.0_dp); o = o + nze
    pg_c = max(reshape(v(o + 1:o + nze), [n_z, n_en]), 0.0_dp); o = o + nze
    pl_v = max(reshape(v(o + 1:o + nze), [n_z, n_en]), 0.0_dp); o = o + nze
    pg_v = max(reshape(v(o + 1:o + nze), [n_z, n_en]), 0.0_dp); o = o + nze
    pp_lc = max(reshape(v(o + 1:o + nze), [n_z, n_en]), 0.0_dp); o = o + nze
    pp_gc = max(reshape(v(o + 1:o + nze), [n_z, n_en]), 0.0_dp); o = o + nze
    pp_lv = max(reshape(v(o + 1:o + nze), [n_z, n_en]), 0.0_dp); o = o + nze
    pp_gv = max(reshape(v(o + 1:o + nze), [n_z, n_en]), 0.0_dp); o = o + nze
    s_l(1, :) = cmplx(0.0_dp, v(o + 1:o + n_en), dp); o = o + n_en
    s_g(1, :) = cmplx(0.0_dp, v(o + 1:o + n_en), dp)
  end subroutine unpack_x

  !-----------------------------------------------------------------------------
  ! Solve the small dense system a x = b (Gaussian elimination, partial pivot)
  !-----------------------------------------------------------------------------
  subroutine solve_small(m, a_in, b_in, x)
    integer, intent(in) :: m
    real(dp), intent(in) :: a_in(n_and, n_and), b_in(n_and)
    real(dp), intent(out) :: x(n_and)

    real(dp) :: a(n_and, n_and), b(n_and), f, tmp(n_and), tb
    integer :: i1, j1, k1, ip1

    a = a_in
    b = b_in
    x = 0.0_dp
    do k1 = 1, m
      ip1 = k1
      do i1 = k1 + 1, m
        if (abs(a(i1, k1)) > abs(a(ip1, k1))) ip1 = i1
      end do
      if (ip1 /= k1) then
        tmp(:) = a(k1, :)
        a(k1, :) = a(ip1, :)
        a(ip1, :) = tmp(:)
        tb = b(k1)
        b(k1) = b(ip1)
        b(ip1) = tb
      end if
      if (abs(a(k1, k1)) < 1.0e-300_dp) cycle
      do i1 = k1 + 1, m
        f = a(i1, k1) / a(k1, k1)
        a(i1, k1:m) = a(i1, k1:m) - f * a(k1, k1:m)
        b(i1) = b(i1) - f * b(k1)
      end do
    end do
    do i1 = m, 1, -1
      if (abs(a(i1, i1)) < 1.0e-300_dp) cycle
      f = b(i1)
      do j1 = i1 + 1, m
        f = f - a(i1, j1) * x(j1)
      end do
      x(i1) = f / a(i1, i1)
    end do
  end subroutine solve_small

  !-----------------------------------------------------------------------------
  ! Energy-integrated electron (conduction chain) and hole (valence chain)
  ! electric currents at the bonds, mA/cm^2, positive along +z.
  !-----------------------------------------------------------------------------
  subroutine currents(je, jh)
    real(dp), intent(out) :: je(n_z - 1), jh(n_z - 1)
    integer :: bb

    do bb = 1, n_z - 1
      je(bb) = cf * sum(rci(1:bb, :))
      jh(bb) = cf * sum(rvi(1:bb, :))
    end do
  end subroutine currents

  !-----------------------------------------------------------------------------
  ! Spectral current j(E,z) [mA/(cm^2 eV)] at the bonds, from the cumulative sum
  ! of the spectral in-scattering rates.
  !-----------------------------------------------------------------------------
  subroutine write_spectrum(lab)
    character(len=1), intent(in) :: lab

    real(dp) :: accc, accv, cfs
    integer :: ie, bb, fus

    cfs = -q_e / (2.0_dp * pi * hbar_evs) * 1.0e3_dp
    open(newunit=fus, file='hetero_interface_spec_'//lab//'.dat', &
        & status='replace', action='write')
    write(fus, '(a)') '# band(1=c,2=v)  z[nm]  E[eV]  j(E,z) [mA/(cm^2 eV)]'
    do ie = 1, n_en
      accc = 0.0_dp
      accv = 0.0_dp
      do bb = 1, n_z - 1
        accc = accc + rc(bb, ie)
        accv = accv + rv(bb, ie)
        if (e_grid(ie) >= 1.2_dp) then
          write(fus, '(i2,f8.3,f9.4,es14.6)') 1, (real(bb, dp) - 0.5_dp) * dz_nm, &
              & e_grid(ie), cfs * accc
        end if
        if (e_grid(ie) <= 0.4_dp) then
          write(fus, '(i2,f8.3,f9.4,es14.6)') 2, (real(bb, dp) - 0.5_dp) * dz_nm, &
              & e_grid(ie), cfs * accv
        end if
      end do
      write(fus, '()')
    end do
    close(fus)
  end subroutine write_spectrum

end program hetero_interface
