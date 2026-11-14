!===============================================================================
! test bulk_srh_rate
!
! PURPOSE
!   Compare the Shockley-Read-Hall (SRH) volume recombination rate of a bulk
!   semiconductor with a single midgap defect level, computed four ways:
!
!     r_schenk       semiclassical rate, multiphonon capture coefficients,
!                    Fermi statistics (reference), approach in paper
!                    [Schenk1992Improved]
!     r_oneshot_bare NEGF, one-shot (Born approximation) with NEITHER Green's
!                    function renormalised: the defect spectral function stays
!                    a delta function at eps_d and the band Green's functions
!                    are the bare ones (contacts only).
!     r_oneshot_band NEGF, one iteration of the self-consistency loop. The
!                    defect Green's function is NOT renormalised (delta
!                    function at eps_d). The band Green's functions ARE
!                    renormalised, once, by the defect self-energy.
!     r_full         NEGF, self-consistent until convergence. BOTH Green's
!                functions are renormalised: the defect one by the Dyson
!                equation (its spectral function is broadened by the coupling
!                to the bands), and the band ones by the defect self-energy and
!                by intraband LO-phonon scattering, which provides the cascade
!                capture into the defect.
!
!   The bulk is GaAs-like, the defect density is 1e14 cm^-3,
!   and the bias V is the quasi-Fermi-level splitting mu_c - mu_v = qV.
!   GaAs bulk parameters is taken from Table I of [Aeberhard2019Nonequilibrium].
!   This test is intended for reproducing Fig. 2b of this paper.
!
! PHYSICS OF THE NEGF MODEL
!   Defect level eps_d in the gap couples to the conduction band (energy
!   E >= Ec) and to the valence band (E <= Ev) by multiphonon emission and
!   capture with Huang-Rhys factor S and phonon energy hw.
!   The defect occupation f_d follows from the steady-state balance of the
!   electron and hole capture/emission rates (solve_defect_occupation for a
!   delta-function defect, solve_defect_dyson_keldysh for the renormalised one).
!   The steady-state recombination rate per unit volume is
!       r = rho_d * R_net / hbar.
!
! OUTPUT
!   bulk_srh_rates.dat : columns  V[V]  r_schenk  r_oneshot_bare
!                        r_oneshot_band  r_full  n_iter
!                        rates in m^-3 s^-1, n_iter = SCBA iterations for r_full.
!
! PARAMETERS
!   From Table I of the reference: m_c = 0.067 m0, m_v = 0.1 m0, Ec = 1.42 eV,
!   Ev = 0, eps_d = 0.71 eV, sigma = 1e-15 cm^2, hw = 36 meV, rho_d = 1e14 cm^-3.
!   Not given there (assumed here): T = 300 K, S = 3.5, contact broadening
!   Gamma = 1 meV, LO-phonon strength giving Gamma_ph = 4 meV, and for the
!   renormalised defect: numerical broadening eta_d = 1 meV, energy window
!   +-50 meV around eps_d with 41 points.
!
! CALIBRATION
!   The unknown prefactor of the multiphonon matrix element is fixed so that the
!   thermal capture coefficient equals v_th * sigma, with v_th = sqrt(3kT/m*).
!   The same calibration is used for the semiclassical and the NEGF rates.
!===============================================================================
program bulk_srh_rates

  use ln_defect_srh, only : solve_defect_occupation, &
      & solve_defect_dyson_keldysh, compute_extended_defect_selfenergy
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
  real(dp), parameter :: m_v = 0.1_dp                    ! hole mass / m0
  real(dp), parameter :: ec = 1.42_dp                    ! conduction band edge, eV
  real(dp), parameter :: ev = 0.0_dp                     ! valence band edge, eV
  real(dp), parameter :: eps_d = 0.71_dp                 ! defect level, eV
  real(dp), parameter :: sigma_cm2 = 1.0e-15_dp          ! capture cross section, cm^2
  real(dp), parameter :: hw = 0.036_dp                   ! phonon energy, eV
  real(dp), parameter :: s_hr = 3.5_dp                   ! Huang-Rhys factor (assumed)
  real(dp), parameter :: rho_d = 1.0e14_dp               ! defect density, cm^-3

  !-----------------------------------------------------------------------------
  ! Numerical parameters
  !-----------------------------------------------------------------------------
  integer,  parameter :: n_modes = 40                    ! multiphonon ladder length
  integer,  parameter :: n_en = 1001                     ! energy grid points
  real(dp), parameter :: de = 0.003_dp                   ! energy step, eV (hw/de = 12)
  integer,  parameter :: k_hw = 12                       ! hw in grid points
  integer,  parameter :: n_st = 3200                     ! band states per band
  real(dp), parameter :: dx = 0.25e-3_dp                 ! band-state spacing, eV
  real(dp), parameter :: gam = 1.0e-3_dp                 ! contact broadening, eV (assumed)
  real(dp), parameter :: gam_lo = 4.0e-3_dp              ! LO broadening target, eV (assumed)
  integer,  parameter :: max_it = 200                    ! SCBA iteration limit
  integer,  parameter :: min_it = 3                      ! minimum SCBA iterations
  real(dp), parameter :: tol_it = 1.0e-9_dp              ! relative tolerance on rate
  real(dp), parameter :: mix = 0.7_dp                    ! mixing of phonon self-energy
  integer,  parameter :: n_omega = 41                    ! defect window points (assumed)
  real(dp), parameter :: omega_hw = 0.05_dp              ! defect window half-width, eV (assumed)
  real(dp), parameter :: eta_d = 1.0e-3_dp               ! defect broadening, eV (assumed)

  !-----------------------------------------------------------------------------
  ! Bias scan
  !-----------------------------------------------------------------------------
  integer,  parameter :: n_volt = 10
  real(dp), parameter :: v0 = 0.6_dp                     ! first bias, V
  real(dp), parameter :: dvolt = 0.05_dp                 ! bias step, V
  real(dp), parameter :: mu_mid = 0.71_dp                ! midpoint of mu_c, mu_v, eV

  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  real(dp) :: kbt, pref_c, pref_v, neff_c, neff_v, vth_n, vth_p
  real(dp) :: cbar_n, cbar_p, c0n, c0p, avg_n, avg_p, e_l, nd_l
  real(dp) :: m0_dc, m0_dv, nb, d2_c, d2_v, ni
  real(dp) :: e_grid(n_en), w_quad(n_en)
  real(dp) :: eps_cs(n_st), w_cs(n_st), eps_vs(n_st), w_vs(n_st)
  real(dp) :: volt, mu_c, mu_v, af, bal
  real(dp) :: r_sc, r_bare, r_band, r_full, r_def_sc, f_sc
  real(dp) :: res(6, 0:n_volt)
  integer :: n_it, n_it_one
  type(TDefectSpec) :: spec
  type(TDefectEnsemble) :: ens_one, ens_full
  integer :: i, l, iv, fu

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

  !=============================================================================
  ! STEP 2: Calibrate the multiphonon prefactor (Schenk semiclassical)
  !   <c_n> = v_th sigma   and   <c_p> = v_th sigma
  !=============================================================================
  avg_n = 0.0_dp
  avg_p = 0.0_dp
  do l = 1, n_modes
    e_l = eps_d + real(l, dp) * hw
    if (e_l >= ec) then
      nd_l = pref_c * sqrt(e_l - ec)
      avg_n = avg_n + shape_n(l) * nd_l * exp(-(e_l - ec) / kbt)
    end if
    e_l = eps_d - real(l, dp) * hw
    if (e_l <= ev) then
      nd_l = pref_v * sqrt(ev - e_l)
      avg_p = avg_p + shape_p(l) * nd_l * exp(-(ev - e_l) / kbt)
    end if
  end do
  c0n = cbar_n / (avg_n / neff_c)
  c0p = cbar_p / (avg_p / neff_v)

  ! Same calibration, expressed through the libnegf kernel (NEGF side)
  call init_defect_spec(spec, 1, [1], [1.0_dp], n_modes)
  call configure_spec(spec, 1.0_dp, 1.0_dp, 1.0_dp)
  avg_n = 0.0_dp
  avg_p = 0.0_dp
  do l = 1, n_modes
    e_l = eps_d + spec%e_shift(l)
    if (e_l >= ec) then
      nd_l = pref_c * sqrt(e_l - ec)
      avg_n = avg_n + 2.0_dp * pi / hbar_evs * spec%m_dc_capt(l) * nd_l &
          & * exp(-(e_l - ec) / kbt)
    end if
    e_l = eps_d - spec%e_shift(l)
    if (e_l <= ev) then
      nd_l = pref_v * sqrt(ev - e_l)
      avg_p = avg_p + 2.0_dp * pi / hbar_evs * spec%m_dv_em(l) * nd_l &
          & * exp(-(ev - e_l) / kbt)
    end if
  end do
  avg_n = avg_n / neff_c
  avg_p = avg_p / neff_v
  m0_dc = cbar_n / avg_n
  m0_dv = cbar_p / avg_p
  call configure_spec(spec, m0_dc, m0_dv, 1.0_dp)

  ! LO-phonon coupling of the bands (Einstein mode, see header)
  nb = 1.0_dp / (exp(hw / kbt) - 1.0_dp)
  d2_c = gam_lo / ((nb + 1.0_dp) * 2.0_dp * pi * pref_c * sqrt(hw))
  d2_v = gam_lo / ((nb + 1.0_dp) * 2.0_dp * pi * pref_v * sqrt(hw))

  ni = sqrt(neff_c * neff_v) * exp(-(ec - ev) / (2.0_dp * kbt))
  write(*, '(a)') '--- parameters ------------------------------------------'
  write(*, '(a,es12.4,a,es12.4)') ' N_c,eff [cm^-3] = ', neff_c, &
      & '   N_v,eff = ', neff_v
  write(*, '(a,es12.4)')          ' n_i [cm^-3]     = ', ni
  write(*, '(a,es12.4,a,es12.4)') ' v_th*sigma (n,p) [cm^3/s] = ', cbar_n, &
      & '   ', cbar_p
  write(*, '(a,es12.4,a,es12.4)') ' M0_dc = ', m0_dc, '   M0_dv = ', m0_dv
  write(*, '(a,es12.4,a,es12.4,a,f8.5)') ' D2_c = ', d2_c, '   D2_v = ', d2_v, &
      & '   n_B = ', nb

  !=============================================================================
  ! STEP 3: Energy grid, band states and defect ensembles
  !=============================================================================
  do i = 1, n_en
    e_grid(i) = eps_d + de * real(i - (n_en + 1) / 2, dp)
    w_quad(i) = de
  end do

  do i = 1, n_st
    af = (real(i, dp) - 0.5_dp) * dx
    eps_cs(i) = ec + af
    w_cs(i) = pref_c * sqrt(af) * dx
    eps_vs(i) = ev - af
    w_vs(i) = pref_v * sqrt(af) * dx
  end do

  ! One-shot: defect Green's function not renormalised (delta function)
  call init_defect_ensemble(ens_one, 1, 1)
  call init_defect_spec(ens_one%defects(1), 1, [1], [1.0_dp], n_modes)
  call configure_spec(ens_one%defects(1), m0_dc, m0_dv, rho_d)

  ! Full: defect Green's function renormalised (Dyson-Keldysh)
  call init_defect_ensemble(ens_full, 1, 1)
  call init_defect_spec(ens_full%defects(1), 1, [1], [1.0_dp], n_modes)
  call configure_spec(ens_full%defects(1), m0_dc, m0_dv, rho_d)
  ens_full%defects(1)%mode = DEFECT_MODE_DYSON
  ens_full%defects(1)%n_omega = n_omega
  ens_full%defects(1)%omega_halfwidth = omega_hw
  ens_full%defects(1)%eta_d = eta_d

  !=============================================================================
  ! STEP 4: Bias scan, all three methods
  !=============================================================================
  do iv = 0, n_volt
    volt = v0 + real(iv, dp) * dvolt
    mu_c = mu_mid + 0.5_dp * volt
    mu_v = mu_mid - 0.5_dp * volt

    ! ---- (a) semiclassical Schenk rate --------------------------------------
    call schenk_fermi(mu_c, mu_v, r_def_sc, f_sc)
    r_sc = 1.0e6_dp * rho_d * r_def_sc

    ! ---- (b) NEGF one-shot, bare: delta-function defect, bare bands (one
    !          pass, nothing renormalised), no LO phonons ----------------------
    call negf_scba(mu_c, mu_v, ens_one, .false., .false., 1, r_bare, n_it_one)
    bal = abs(ens_one%defects(1)%rate_net_el - ens_one%defects(1)%rate_net_hl) &
        & / abs(ens_one%defects(1)%rate_net_el)
    if (bal > 1.0e-9_dp) then
      write(*, '(a,f6.2,a,es12.4)') 'FAIL: one-shot bare el/hl balance at V=', &
          & volt, ' rel err=', bal
      error stop 1
    end if

    ! ---- (c) NEGF one-shot, band: delta-function defect, bands renormalised
    !          once by the defect self-energy, no LO phonons --------------------
    call negf_scba(mu_c, mu_v, ens_one, .false., .false., 2, r_band, n_it_one)
    bal = abs(ens_one%defects(1)%rate_net_el - ens_one%defects(1)%rate_net_hl) &
        & / abs(ens_one%defects(1)%rate_net_el)
    if (bal > 1.0e-9_dp) then
      write(*, '(a,f6.2,a,es12.4)') 'FAIL: one-shot band el/hl balance at V=', &
          & volt, ' rel err=', bal
      error stop 1
    end if

    ! ---- (d) NEGF full: self-consistent, renormalised defect (Dyson) and
    !          bands (defect self-energy + LO phonons) -------------------------
    call negf_scba(mu_c, mu_v, ens_full, .true., .true., max_it, r_full, n_it)
    if (n_it >= max_it) then
      write(*, '(a,f6.2)') 'FAIL: full SCBA did not converge at V=', volt
      error stop 1
    end if

    res(:, iv) = [volt, r_sc, r_bare, r_band, r_full, real(n_it, dp)]
  end do

  call destruct_defect_spec(spec)
  call destruct_defect_ensemble(ens_one)
  call destruct_defect_ensemble(ens_full)

  !=============================================================================
  ! STEP 5: Output
  !=============================================================================
  open(newunit=fu, file='bulk_srh_rates.dat', status='replace', action='write')
  write(fu, '(a)') '# SRH volume recombination rate for bulk semiconductor with single midgap defect'
  write(fu, '(a)') '# rates in m^-3 s^-1'
  write(fu, '(a)') '# Voltage[V]     r_schenk        r_oneshot_bare        r_oneshot_band        r_full      n_iter'
  do iv = 0, n_volt
    write(fu, '(f9.3,4es18.8e3,f6.0)') res(:, iv)
  end do
  close(fu)

  write(*, '(a)') '--- results ----------------------------------------------'
  write(*, '(a)') '   V[V]     r_schenk        r_oneshot_bare        r_oneshot_band        r_full      n_iter'
  do iv = 0, n_volt
    write(*, '(f9.3,4es18.8e3,f6.0)') res(:, iv)
  end do

  !=============================================================================
  ! STEP 6: Checks
  !=============================================================================
  do iv = 1, n_volt
    if (res(4, iv) <= res(4, iv - 1)) then
      write(*, '(a,f6.2)') 'FAIL: full rate not monotonic in V at V=', res(1, iv)
      error stop 1
    end if
  end do

  write(*, '(a)') 'test f90bulk_srh_rate finished, data in bulk_srh_rates.dat'

contains

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
  ! 3D parabolic density of states N(E) = pref * sqrt(E), [1/(eV cm^3)],
  ! E measured from the band edge in eV
  !-----------------------------------------------------------------------------
  pure function dos_pref(m_rel) result(pref)
    real(dp), intent(in) :: m_rel
    real(dp) :: pref
    pref = (2.0_dp * m_rel * m_e * q_e / hbar_js**2)**1.5_dp &
        & / (2.0_dp * pi**2) * 1.0e-6_dp
  end function dos_pref

  !-----------------------------------------------------------------------------
  ! Multiphonon line-shape weights for capture into the conduction band (n)
  ! and valence band (p), Schenk Eq. (12) with r_F = 0
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
  ! Semiclassical SRH rate per defect with Fermi statistics (Schenk Eqs. 7-9):
  !   R = (a_n a_p - b_n b_p) / (a_n + b_n + a_p + b_p),  f_d = (a_n + b_p)/sum
  !-----------------------------------------------------------------------------
  subroutine schenk_fermi(muc, muv, r_def, fd)
    real(dp), intent(in) :: muc, muv
    real(dp), intent(out) :: r_def, fd

    real(dp) :: a_n, b_n, a_p, b_p, el, cl, nl, fl, den
    integer :: ll

    a_n = 0.0_dp
    b_n = 0.0_dp
    a_p = 0.0_dp
    b_p = 0.0_dp
    do ll = 1, n_modes
      el = eps_d + real(ll, dp) * hw
      if (el >= ec) then
        cl = c0n * shape_n(ll)
        nl = pref_c * sqrt(el - ec)
        fl = fermi((el - muc) / kbt)
        a_n = a_n + cl * nl * fl
        b_n = b_n + cl * nl * exp(-(el - eps_d) / kbt) * fermi(-(el - muc) / kbt)
      end if
      el = eps_d - real(ll, dp) * hw
      if (el <= ev) then
        cl = c0p * shape_p(ll)
        nl = pref_v * sqrt(ev - el)
        fl = fermi((el - muv) / kbt)
        a_p = a_p + cl * nl * fermi(-(el - muv) / kbt)
        b_p = b_p + cl * nl * exp(-(eps_d - el) / kbt) * fl
      end if
    end do
    den = a_n + b_n + a_p + b_p
    fd = (a_n + b_p) / den
    r_def = (a_n * a_p - b_n * b_p) / den
  end subroutine schenk_fermi

  !-----------------------------------------------------------------------------
  ! Set the defect parameters (Table I) and build the libnegf MPR kernel
  !-----------------------------------------------------------------------------
  subroutine configure_spec(sp, m0c, m0v, dens)
    type(TDefectSpec), intent(inout) :: sp
    real(dp), intent(in) :: m0c, m0v, dens
    sp%eps_d = eps_d
    sp%kbt = kbt
    sp%density_weight = dens
    sp%hbar_omega0 = hw
    sp%huang_rhys_n = s_hr
    sp%huang_rhys_p = s_hr
    sp%rf_rph_ratio_n = 0.0_dp
    sp%rf_rph_ratio_p = 0.0_dp
    sp%m0_dc = m0c
    sp%m0_dv = m0v
    sp%ec_cut = ec
    sp%ev_cut = ev
    sp%mode = DEFECT_MODE_DELTA
    sp%pl_index = 1
    call build_mpr_kernel(sp)
  end subroutine configure_spec

  !-----------------------------------------------------------------------------
  ! NEGF self-consistent (SCBA) cycle for the bulk with a defect ensemble.
  !   - bands: discrete states with contact broadening gam, dressed by the
  !     defect self-energy (libnegf) and, if use_lo, by intraband LO-phonon
  !     self-energies
  !   - defect: occupation and net rate from libnegf; if use_dyson the defect
  !     Green's function is renormalised (Dyson-Keldysh), otherwise it stays a
  !     delta function (the mode of ens_in must be set accordingly)
  ! n_pass passes are made (a pass = solve defect with the current band G, then
  ! update the band self-energies); the cycle stops earlier if the rate has
  ! converged. The rate is evaluated at the beginning of each pass, so n_pass = 2
  ! means: bare bands -> self-energies -> bands renormalised once -> rate.
  ! On return: r_conv = recombination rate [m^-3 s^-1], nit = passes used.
  !-----------------------------------------------------------------------------
  subroutine negf_scba(muc, muv, ens_in, use_dyson, use_lo, n_pass, r_conv, nit)
    real(dp), intent(in) :: muc, muv
    type(TDefectEnsemble), intent(inout) :: ens_in
    logical, intent(in) :: use_dyson, use_lo
    integer, intent(in) :: n_pass
    real(dp), intent(out) :: r_conv
    integer, intent(out) :: nit

    complex(dp) :: gl(1, n_en), gg(1, n_en)
    complex(dp) :: s_l(1, n_en), s_g(1, n_en), s_r(1, n_en)
    real(dp) :: tlc(n_en), tgc(n_en), tlv(n_en), tgv(n_en)
    real(dp) :: pl_c(n_en), pg_c(n_en), pl_v(n_en), pg_v(n_en)
    real(dp) :: new_lc, new_gc, new_lv, new_gv
    real(dp) :: ee, fc, fv, dd, re_c, re_v, wid_c, wid_v
    real(dp) :: sl_c, sg_c, sl_v, sg_v
    real(dp) :: sum_c, sum_v, r_old, r_new
    integer :: ie, js, it, ip, im

    s_l(:, :) = (0.0_dp, 0.0_dp)
    s_g(:, :) = (0.0_dp, 0.0_dp)
    s_r(:, :) = (0.0_dp, 0.0_dp)
    pl_c(:) = 0.0_dp
    pg_c(:) = 0.0_dp
    pl_v(:) = 0.0_dp
    pg_v(:) = 0.0_dp
    r_old = -1.0_dp
    r_new = 0.0_dp
    nit = n_pass

    do it = 1, n_pass

      ! -- dress the band states and sum over them -------------------------
      do ie = 1, n_en
        ee = e_grid(ie)
        fc = fermi((ee - muc) / kbt)
        fv = fermi((ee - muv) / kbt)

        re_c = 0.0_dp
        wid_c = 0.5_dp * gam
        sl_c = pl_c(ie)
        sg_c = pg_c(ie)
        if (ee >= ec) then
          re_c = real(s_r(1, ie), dp)
          wid_c = wid_c - aimag(s_r(1, ie))
          sl_c = sl_c + aimag(s_l(1, ie))
          sg_c = sg_c - aimag(s_g(1, ie))
        end if
        wid_c = wid_c + 0.5_dp * (pl_c(ie) + pg_c(ie))

        re_v = 0.0_dp
        wid_v = 0.5_dp * gam
        sl_v = pl_v(ie)
        sg_v = pg_v(ie)
        if (ee <= ev) then
          re_v = real(s_r(1, ie), dp)
          wid_v = wid_v - aimag(s_r(1, ie))
          sl_v = sl_v + aimag(s_l(1, ie))
          sg_v = sg_v - aimag(s_g(1, ie))
        end if
        wid_v = wid_v + 0.5_dp * (pl_v(ie) + pg_v(ie))

        sum_c = 0.0_dp
        sum_v = 0.0_dp
        do js = 1, n_st
          dd = ee - eps_cs(js) - re_c
          sum_c = sum_c + w_cs(js) / (dd * dd + wid_c * wid_c)
          dd = ee - eps_vs(js) - re_v
          sum_v = sum_v + w_vs(js) / (dd * dd + wid_v * wid_v)
        end do

        tlc(ie) = sum_c * (gam * fc + sl_c)
        tgc(ie) = sum_c * (gam * (1.0_dp - fc) + sg_c)
        tlv(ie) = sum_v * (gam * fv + sl_v)
        tgv(ie) = sum_v * (gam * (1.0_dp - fv) + sg_v)

        gl(1, ie) = cmplx(0.0_dp, tlc(ie) + tlv(ie), dp)
        gg(1, ie) = cmplx(0.0_dp, -(tgc(ie) + tgv(ie)), dp)
      end do

      ! -- defect occupation, net rate and defect self-energy (libnegf) -----
      if (use_dyson) then
        call solve_defect_dyson_keldysh(ens_in%defects(1), e_grid, w_quad, gl, gg)
      else
        call solve_defect_occupation(ens_in%defects(1), e_grid, w_quad, gl, gg)
      end if
      r_new = 1.0e6_dp * ens_in%defects(1)%density_weight &
          & * ens_in%defects(1)%rate_net_el / hbar_evs

      s_l(:, :) = (0.0_dp, 0.0_dp)
      s_g(:, :) = (0.0_dp, 0.0_dp)
      s_r(:, :) = (0.0_dp, 0.0_dp)
      call compute_extended_defect_selfenergy(ens_in, e_grid, w_quad, s_l, s_g, s_r)

      ! -- LO-phonon self-energies from the band-resolved G (mixed) ---------
      if (use_lo) then
        do ie = 1, n_en
          ip = ie + k_hw
          im = ie - k_hw
          new_lc = 0.0_dp
          new_gc = 0.0_dp
          new_lv = 0.0_dp
          new_gv = 0.0_dp
          if (ip <= n_en) then
            new_lc = new_lc + (nb + 1.0_dp) * tlc(ip)
            new_gc = new_gc + nb * tgc(ip)
            new_lv = new_lv + (nb + 1.0_dp) * tlv(ip)
            new_gv = new_gv + nb * tgv(ip)
          end if
          if (im >= 1) then
            new_lc = new_lc + nb * tlc(im)
            new_gc = new_gc + (nb + 1.0_dp) * tgc(im)
            new_lv = new_lv + nb * tlv(im)
            new_gv = new_gv + (nb + 1.0_dp) * tgv(im)
          end if
          pl_c(ie) = (1.0_dp - mix) * pl_c(ie) + mix * d2_c * new_lc
          pg_c(ie) = (1.0_dp - mix) * pg_c(ie) + mix * d2_c * new_gc
          pl_v(ie) = (1.0_dp - mix) * pl_v(ie) + mix * d2_v * new_lv
          pg_v(ie) = (1.0_dp - mix) * pg_v(ie) + mix * d2_v * new_gv
        end do
      end if

      ! -- convergence on the recombination rate ----------------------------
      if (it >= min_it .and. r_old > 0.0_dp) then
        if (abs(r_new - r_old) <= tol_it * abs(r_new)) then
          nit = it
          exit
        end if
      end if
      r_old = r_new
    end do

    r_conv = r_new

  end subroutine negf_scba

end program bulk_srh_rates
