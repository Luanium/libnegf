!===============================================================================
! test pin_defect_current
!
! PURPOSE
!   Recombination current density J(V) of a 100 nm thick GaAs p-i-n diode with
!   a single plane of midgap defects in the middle of the intrinsic region,
!   in the dark, for forward bias V = 0.6 ... 1.1 V. Computed for two doping
!   densities of the contact regions, 1e18 and 1e17 cm^-3, with
!
!     J_semiclassical  standard SRH recombination at the defect position as in
!                drift-diffusion simulators (AFORS-HET): Boltzmann statistics,
!                band edges from the Poisson band profile, capture time
!                constants from v_th * sigma with v_th = sqrt(3kT/m*) (same
!                as in the NEGF calibration), no tunnelling into the gap.
!                Computed for the 1e18 profile.
!     J_ideality2  reference line J0 [exp(V/2kT) - 1] of ideality factor 2,
!                with J0 fixed so that it passes through J_full(1e17) at the
!                highest bias (lowest field, least tunnelling).
!     J_full     NEGF, self-consistent until convergence, with BOTH Green's
!                functions renormalised: the defect one by the Dyson equation
!                (spectral function broadened by the coupling to the bands), and
!                the band ones by the defect self-energy and by intraband
!                LO-phonon scattering (cascade capture). The band Green's
!                functions are resolved in space (1D real-space chain along the
!                growth direction z, 2D transverse momentum), so that the
!                field-induced band tails (Franz-Keldysh) in the gap at the
!                defect position are included, and with them the tunnelling-
!                assisted capture. Reference: [Aeberhard2019Nonequilibrium].
!
!   The expected qualitative behaviour: with the large built-in field of the
!   1e18 doping the NEGF current exceeds the semiclassical one and its slope
!   in V deviates from the semiclassical ideality factor of 2 because of the
!   tunnelling into field-induced band tails; for the 1e17 doping the lower
!   field brings the NEGF characteristic back towards the semiclassical one.
!
!   The results are to be compared with Fig. 3b in the above referenced paper.
!   Because many parameters are not explicitly mentioned in that paper,
!   discrepancies are expected.
!
! DEVICE AND BAND PROFILE
!   100 nm, 1 nm cells: 20 nm p-doped, 60 nm intrinsic, 20 nm n-doped, defect
!   plane in cell 50 (z ~ 50 nm). The electrostatic potential energy U(z) is
!   obtained from a 1D Poisson equation with Boltzmann carrier statistics and
!   flat quasi-Fermi levels, mu_c - mu_v = qV (mu_v = 0). Band edges:
!   Ec(z) = Ec0 + U(z), Ev(z) = Ev0 + U(z), defect level eps_d(z) = eps_d0 + U(z).
!   Conduction and valence band are two decoupled effective-mass chains, each
!   coupled at every site to its own reservoir (broadening gam, quasi-Fermi
!   level mu_c resp. mu_v), as in the bulk case.
!
! NEGF MODEL OF THE DEFECT
!   Level eps_d couples to the conduction band (E > eps_d) and to the valence
!   band (E < eps_d) by multiphonon emission and capture with Huang-Rhys factor
!   S and phonon energy hw. The defect occupation follows from the steady-state
!   balance of electron and hole capture/emission (solve_defect_dyson_keldysh),
!   the volume recombination rate at the defect cell is
!       r = rho_d * R_net / hbar,
!   and the areal recombination current is J = q * r * dz (the defect volume
!   density acts over one cell of thickness dz).
!
! OUTPUT
!   pin_defect_current.dat : columns  V[V]  J_semiclassical  J_full(1e18)
!                            J_full(1e17)  J_ideality2  n_iter(1e18)  n_iter(1e17)
!                            J in mA/cm^2, n_iter = SCBA iterations.
!                            J_semiclassical = semiclassical (drift-diffusion-like)
!                            local SRH current, J_ideality2 = J0 [exp(V/2kT) - 1].
!
! PARAMETERS
!   From the reference: m_c = 0.067 m0, m_v = 0.22 m0, eps_d = 0.75 eV,
!   rho_d = 1e14 cm^-3, sigma_n = sigma_p = 1e-14 cm^2, hw = 36 meV, Ec = 1.42 eV,
!   Ev = 0. Derived from these (not assigned): effective DOS Nc, Nv and thermal
!   velocities v_th = sqrt(3kT/m*). Not given there (assumed here): T = 300 K,
!   S = 3.5, layer widths (20 nm p, 60 nm intrinsic, 20 nm n), relative
!   permittivity 12.9, contact broadening gam = 1 meV, LO-phonon strength giving
!   Gamma_ph = 4 meV, 1 nm discretisation and 1 nm defect layer, numerical
!   broadening eta_d = 1 meV and window +-50 meV (41 points) for the
!   renormalised defect, and the depth tun_win = 50 meV of the field-induced
!   tail states below the local band edge that are coupled to the defect.
!
! CALIBRATION
!   The unknown prefactor of the multiphonon matrix element is fixed (flat-band
!   bulk) so that the thermal capture coefficient equals v_th * sigma, with
!   v_th = sqrt(3kT/m*). For the NEGF side the band DOS entering the average is
!   the one actually sampled by the discrete chain with the transverse
!   quadrature (chain_dos), so that the NEGF rate in flat-band, quasi-
!   equilibrium conditions equals the SRH rate. The same v_th*sigma, and Nc, Nv
!   from the effective masses, are used in the semiclassical rate.
!
! KNOWN DIFFERENCE TO THE REFERENCE
!   J_full(1e18), J_full(1e17) and the ideality-2 line reproduce Fig. 3b well.
!   J_semiclassical has the correct slope (ideality factor 2) but is a
!   factor ~2 above the drift-diffusion curve of the paper, which was computed
!   with AFORS-HET and its internal, unpublished parameters (v_th, Nc, Nv,
!   defect layer thickness). No parameter was tuned to remove this offset, in
!   order to keep the semiclassical and NEGF rates mutually consistent.
!===============================================================================
program pin_defect_current

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
  real(dp), parameter :: eps0 = 8.8541878128e-12_dp      ! F/m

  !-----------------------------------------------------------------------------
  ! Material and defect parameters
  !-----------------------------------------------------------------------------
  real(dp), parameter :: temp = 300.0_dp                 ! K (assumed)
  real(dp), parameter :: m_c = 0.067_dp                  ! electron mass / m0
  real(dp), parameter :: m_v = 0.22_dp                   ! hole mass / m0
  real(dp), parameter :: ec0 = 1.42_dp                   ! conduction band edge, eV
  real(dp), parameter :: ev0 = 0.0_dp                    ! valence band edge, eV
  real(dp), parameter :: eps_d0 = 0.75_dp                ! defect level (flat band), eV
  real(dp), parameter :: sigma_cm2 = 1.0e-14_dp          ! capture cross section, cm^2
  real(dp), parameter :: hw = 0.036_dp                   ! phonon energy, eV
  real(dp), parameter :: s_hr = 3.5_dp                   ! Huang-Rhys factor (assumed)
  real(dp), parameter :: rho_d = 1.0e14_dp               ! defect density, cm^-3
  real(dp), parameter :: eps_r = 12.9_dp                 ! relative permittivity (assumed)

  !-----------------------------------------------------------------------------
  ! Device
  !-----------------------------------------------------------------------------
  integer,  parameter :: n_z = 100                       ! cells
  real(dp), parameter :: dz_nm = 1.0_dp                  ! cell thickness, nm (assumed)
  integer,  parameter :: nz_p = 20                       ! p-doped cells (assumed)
  integer,  parameter :: nz_n = 20                       ! n-doped cells (assumed)
  integer,  parameter :: i_def = 50                      ! defect cell
  real(dp), parameter :: ndop(2) = [1.0e18_dp, 1.0e17_dp] ! doping, cm^-3

  !-----------------------------------------------------------------------------
  ! Numerical parameters
  !-----------------------------------------------------------------------------
  integer,  parameter :: n_modes = 40                    ! multiphonon ladder length
  integer,  parameter :: n_en = 1041                     ! energy grid points
  real(dp), parameter :: de = 0.003_dp                   ! energy step, eV (hw/de = 12)
  integer,  parameter :: k_hw = 12                       ! hw in grid points
  integer,  parameter :: n_kt = 100                      ! transverse energy nodes
  real(dp), parameter :: gam = 1.0e-3_dp                 ! contact broadening, eV (assumed)
  real(dp), parameter :: gam_lo = 4.0e-3_dp              ! LO broadening target, eV (assumed)
  integer,  parameter :: max_it = 200                    ! SCBA iteration limit
  integer,  parameter :: min_it = 3                      ! minimum SCBA iterations
  real(dp), parameter :: tol_it = 1.0e-6_dp              ! relative tolerance on rate
  real(dp), parameter :: mix = 0.7_dp                    ! mixing of phonon self-energy
  integer,  parameter :: n_omega = 41                    ! defect window points (assumed)
  real(dp), parameter :: omega_hw = 0.05_dp              ! defect window half-width, eV (assumed)
  real(dp), parameter :: eta_d = 1.0e-3_dp               ! defect broadening, eV (assumed)
  real(dp), parameter :: tun_win = 0.05_dp               ! max depth of tail states below the local band edge coupled to the defect, eV (assumed)

  !-----------------------------------------------------------------------------
  ! Bias scan
  !-----------------------------------------------------------------------------
  integer,  parameter :: n_volt = 5
  real(dp), parameter :: v0 = 0.6_dp                     ! first bias, V
  real(dp), parameter :: dvolt = 0.1_dp                 ! bias step, V

  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  real(dp) :: kbt, pref_c, pref_v, neff_c, neff_v, vth_n, vth_p
  real(dp) :: cbar_n, cbar_p, c0n, c0p, avg_n, avg_p, e_l, nd_l
  real(dp) :: m0_dc, m0_dv, nb, d2_c, d2_v
  real(dp) :: dz_cm, t_c, t_v, pk_c, pk_v
  real(dp) :: u(n_z), eps_dl, ec_lo, ev_hi
  real(dp) :: volt, mu_c, mu_v, r_def_sc, r_full
  real(dp) :: j_sc(0:n_volt), j_full(2, 0:n_volt), nit_arr(2, 0:n_volt)
  real(dp) :: j_ref(0:n_volt), j0_ref
  real(dp), allocatable :: e_grid(:), w_quad(:)
  real(dp), allocatable :: tlc(:, :), tgc(:, :), tlv(:, :), tgv(:, :)
  real(dp), allocatable :: pl_c(:, :), pg_c(:, :), pl_v(:, :), pg_v(:, :)
  complex(dp), allocatable :: s_l(:, :), s_g(:, :), s_r(:, :)
  integer :: n_it
  type(TDefectSpec) :: spec
  type(TDefectEnsemble) :: ens
  integer :: i, l, iv, id, fu

  !=============================================================================
  ! STEP 1: Constants and band quantities
  !=============================================================================
  kbt = kb_ev * temp
  pref_c = dos_pref(m_c)
  pref_v = dos_pref(m_v)
  neff_c = pref_c * 0.5_dp * sqrt(pi) * kbt**1.5_dp      ! effective DOS, cm^-3
  neff_v = pref_v * 0.5_dp * sqrt(pi) * kbt**1.5_dp
  ! thermal velocities, shared by NEGF calibration and semiclassical SRH
  vth_n = 100.0_dp * sqrt(3.0_dp * kbt * q_e / (m_c * m_e))   ! cm/s
  vth_p = 100.0_dp * sqrt(3.0_dp * kbt * q_e / (m_v * m_e))
  cbar_n = vth_n * sigma_cm2                             ! cm^3/s
  cbar_p = vth_p * sigma_cm2

  dz_cm = dz_nm * 1.0e-7_dp
  ! hopping of the effective-mass chains, eV
  t_c = hbar_js**2 / (2.0_dp * m_c * m_e * (dz_nm * 1.0e-9_dp)**2) / q_e
  t_v = hbar_js**2 / (2.0_dp * m_v * m_e * (dz_nm * 1.0e-9_dp)**2) / q_e
  ! transverse (2D, with spin) density of states m/(pi hbar^2), cm^-2 eV^-1
  pk_c = m_c * m_e / (pi * hbar_js**2) * q_e * 1.0e-4_dp
  pk_v = m_v * m_e / (pi * hbar_js**2) * q_e * 1.0e-4_dp

  !=============================================================================
  ! STEP 2: Calibrate the multiphonon prefactor (flat-band bulk)
  !   <c_n> = v_th sigma   and   <c_p> = v_th sigma
  !=============================================================================
  avg_n = 0.0_dp
  avg_p = 0.0_dp
  do l = 1, n_modes
    e_l = eps_d0 + real(l, dp) * hw
    if (e_l >= ec0) then
      nd_l = pref_c * sqrt(e_l - ec0)
      avg_n = avg_n + shape_n(l) * nd_l * exp(-(e_l - ec0) / kbt)
    end if
    e_l = eps_d0 - real(l, dp) * hw
    if (e_l <= ev0) then
      nd_l = pref_v * sqrt(ev0 - e_l)
      avg_p = avg_p + shape_p(l) * nd_l * exp(-(ev0 - e_l) / kbt)
    end if
  end do
  c0n = cbar_n / (avg_n / neff_c)
  c0p = cbar_p / (avg_p / neff_v)

  ! Same calibration, expressed through the libnegf kernel (NEGF side). The DOS
  ! is the one actually sampled by the discrete chain + transverse quadrature
  ! (chain_dos), so that the NEGF flat-band rate equals the SRH one by design.
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

  ! LO-phonon coupling of the bands (Einstein mode, see header)
  nb = 1.0_dp / (exp(hw / kbt) - 1.0_dp)
  d2_c = gam_lo / ((nb + 1.0_dp) * 2.0_dp * pi * pref_c * sqrt(hw))
  d2_v = gam_lo / ((nb + 1.0_dp) * 2.0_dp * pi * pref_v * sqrt(hw))

  write(*, '(a)') '--- parameters ------------------------------------------'
  write(*, '(a,es12.4,a,es12.4)') ' N_c,eff [cm^-3] = ', neff_c, &
      & '   N_v,eff = ', neff_v
  write(*, '(a,es12.4,a,es12.4)') ' v_th*sigma (n,p) [cm^3/s] = ', cbar_n, &
      & '   ', cbar_p
  write(*, '(a,es12.4,a,es12.4)') ' M0_dc = ', m0_dc, '   M0_dv = ', m0_dv
  write(*, '(a,es12.4,a,es12.4)') ' t_c [eV] = ', t_c, '   t_v [eV] = ', t_v

  !=============================================================================
  ! STEP 3: Defect ensemble (renormalised defect, Dyson-Keldysh) and arrays
  !=============================================================================
  call init_defect_ensemble(ens, 1, 1)
  call init_defect_spec(ens%defects(1), 1, [1], [1.0_dp], n_modes)
  call configure_spec(ens%defects(1), m0_dc, m0_dv, rho_d)
  ens%defects(1)%mode = DEFECT_MODE_DYSON
  ens%defects(1)%n_omega = n_omega
  ens%defects(1)%omega_halfwidth = omega_hw
  ens%defects(1)%eta_d = eta_d

  allocate(e_grid(n_en), w_quad(n_en))
  allocate(tlc(n_z, n_en), tgc(n_z, n_en), tlv(n_z, n_en), tgv(n_z, n_en))
  allocate(pl_c(n_z, n_en), pg_c(n_z, n_en), pl_v(n_z, n_en), pg_v(n_z, n_en))
  allocate(s_l(1, n_en), s_g(1, n_en), s_r(1, n_en))
  w_quad(:) = de

  !=============================================================================
  ! STEP 4: Bias scan for both dopings
  !=============================================================================
  do id = 1, 2
    do iv = 0, n_volt
      volt = v0 + real(iv, dp) * dvolt
      mu_c = volt
      mu_v = 0.0_dp

      ! band profile from Poisson
      call solve_poisson(volt, ndop(id), u)
      eps_dl = eps_d0 + u(i_def)
      ec_lo = ec0 + minval(u) - 0.1_dp
      ev_hi = ev0 + maxval(u) + 0.1_dp
      do i = 1, n_en
        e_grid(i) = eps_dl + de * real(i - (n_en + 1) / 2, dp)
      end do

      ! ---- semiclassical local SRH rate at the defect (1e18 profile) --------
      if (id == 1) then
        call srh_standard(mu_c, mu_v, ec0 + u(i_def), ev0 + u(i_def), eps_dl, &
            & r_def_sc)
        j_sc(iv) = q_e * r_def_sc * dz_cm * 1.0e3_dp
      end if

      ! ---- NEGF full ---------------------------------------------------------
      ens%defects(1)%eps_d = eps_dl
      ! capture/emission resonances only for band states and Franz-Keldysh
      ! tail states within tun_win of the local band edges at the defect
      ens%defects(1)%ec_cut = ec0 + u(i_def) - tun_win
      ens%defects(1)%ev_cut = ev0 + u(i_def) + tun_win
      call negf_scba(mu_c, mu_v, r_full, n_it)
      if (n_it >= max_it) then
        write(*, '(a,f6.2,a,es10.2)') 'FAIL: full SCBA did not converge at V=', &
            & volt, '  doping=', ndop(id)
        error stop 1
      end if
      j_full(id, iv) = q_e * r_full * dz_cm * 1.0e3_dp
      nit_arr(id, iv) = real(n_it, dp)
      write(*, '(a,es9.2,a,f6.2,a,es14.6,a,i0)') ' doping=', ndop(id), '  V=', &
          & volt, '  J_full=', j_full(id, iv), '  n_iter=', n_it
    end do
  end do

  call destruct_defect_ensemble(ens)

  ! ---- ideality-2 reference line J0 [exp(V/2kT) - 1] ------------------------
  ! J0 is fixed so that the line passes through the NEGF current of the lowest
  ! doping (1e17, smallest built-in field, closest to pure ideality-2 behaviour
  ! without tunnelling) at the highest bias; the semiclassical curve stays
  ! independent.
  j0_ref = j_full(2, n_volt) / (exp((v0 + real(n_volt, dp) * dvolt) / (2.0_dp * kbt)) &
      & - 1.0_dp)
  do iv = 0, n_volt
    j_ref(iv) = j0_ref * (exp((v0 + real(iv, dp) * dvolt) / (2.0_dp * kbt)) &
        & - 1.0_dp)
  end do

  !=============================================================================
  ! STEP 5: Output
  !=============================================================================
  open(newunit=fu, file='pin_defect_current.dat', status='replace', action='write')
  write(fu, '(a)') '# Recombination current of a GaAs p-i-n diode with a midgap defect plane'
  write(fu, '(a)') '# current density in mA/cm^2'
  write(fu, '(a)') '# Voltage[V]  J_semiclassical  J_full(1e18)  J_full(1e17)  J_ideality2  n_iter(1e18)  n_iter(1e17)'
  do iv = 0, n_volt
    write(fu, '(f9.3,4es18.8e3,2f6.0)') v0 + real(iv, dp) * dvolt, j_sc(iv), &
        & j_full(1, iv), j_full(2, iv), j_ref(iv), nit_arr(1, iv), nit_arr(2, iv)
  end do
  close(fu)

  write(*, '(a)') '--- results ----------------------------------------------'
  write(*, '(a)') '   V[V]   J_semiclassical      J_full(1e18)      J_full(1e17)       J_ideality2  n_it(1e18) n_it(1e17)'
  do iv = 0, n_volt
    write(*, '(f9.3,4es18.8e3,2f6.0)') v0 + real(iv, dp) * dvolt, j_sc(iv), &
        & j_full(1, iv), j_full(2, iv), j_ref(iv), nit_arr(1, iv), nit_arr(2, iv)
  end do

  !=============================================================================
  ! STEP 6: Checks
  !=============================================================================
  do id = 1, 2
    do iv = 1, n_volt
      if (j_full(id, iv) <= j_full(id, iv - 1)) then
        write(*, '(a,i0,a,f6.2)') 'FAIL: full current not monotonic in V, doping index ', &
            & id, ' at V=', v0 + real(iv, dp) * dvolt
        error stop 1
      end if
    end do
  end do

  write(*, '(a)') 'test f90pin_defect_current finished, data in pin_defect_current.dat'

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
  ! DOS [1/(eV cm^3)] of the discrete chain with transverse midpoint quadrature
  ! and contact broadening gam, exactly as sampled by band_energy for a
  ! flat-band, uniform region. x_in = energy into the band measured from the
  ! local band edge (eV). ib = 1: conduction, ib = 2: valence.
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
  ! Standard semiclassical SRH recombination rate (as in drift-diffusion
  ! simulators such as AFORS-HET) for a single defect level with Boltzmann
  ! carrier statistics and flat quasi-Fermi levels muc, muv at the defect
  ! position, for local band edges ec_l, ev_l and local defect level ed_l:
  !   R = (n p - n1 p1) / (tau_p (n + n1) + tau_n (p + p1))   [cm^-3 s^-1]
  !   tau_n = 1/(v_th,n sigma_n rho_d),  tau_p = 1/(v_th,p sigma_p rho_d)
  !   n = Nc exp((muc - ec_l)/kT),  p = Nv exp((ev_l - muv)/kT)
  !   n1 = Nc exp((ed_l - ec_l)/kT),  p1 = Nv exp((ev_l - ed_l)/kT)
  !-----------------------------------------------------------------------------
  subroutine srh_standard(muc, muv, ec_l, ev_l, ed_l, r_vol)
    real(dp), intent(in) :: muc, muv, ec_l, ev_l, ed_l
    real(dp), intent(out) :: r_vol

    real(dp) :: n_c, p_v, n1, p1, tau_n, tau_p

    n_c = neff_c * exp((muc - ec_l) / kbt)
    p_v = neff_v * exp((ev_l - muv) / kbt)
    n1 = neff_c * exp((ed_l - ec_l) / kbt)
    p1 = neff_v * exp((ev_l - ed_l) / kbt)
    tau_n = 1.0_dp / (cbar_n * rho_d)
    tau_p = 1.0_dp / (cbar_p * rho_d)
    r_vol = (n_c * p_v - n1 * p1) / (tau_p * (n_c + n1) + tau_n * (p_v + p1))
  end subroutine srh_standard

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
  ! 1D Poisson equation with Boltzmann statistics and flat quasi-Fermi levels
  ! mu_v = 0, mu_c = vb. Returns the electrostatic potential energy U(z) [eV]
  ! (Ec = Ec0 + U, Ev = Ev0 + U). The first and last cell are charge neutral
  ! contact cells.
  !-----------------------------------------------------------------------------
  subroutine solve_poisson(vb, ndop_cm3, uvec)
    real(dp), intent(in) :: vb, ndop_cm3
    real(dp), intent(out) :: uvec(n_z)

    real(dp) :: dopn(n_z), rhs(n_z), bb(n_z), cp(n_z), dq(n_z), du(n_z)
    real(dp) :: nc_m3, nv_m3, mu_cl, mu_vl, dzm, cpois, off, pp, nn, den, step
    integer :: it, ip

    dzm = dz_nm * 1.0e-9_dp
    cpois = q_e / (eps_r * eps0)
    off = 1.0_dp / dzm**2
    nc_m3 = neff_c * 1.0e6_dp
    nv_m3 = neff_v * 1.0e6_dp
    mu_vl = 0.0_dp
    mu_cl = vb

    dopn(:) = 0.0_dp
    dopn(1:nz_p) = -ndop_cm3 * 1.0e6_dp
    dopn(n_z - nz_n + 1:n_z) = ndop_cm3 * 1.0e6_dp

    uvec(1) = mu_vl + kbt * log(ndop_cm3 * 1.0e6_dp / nv_m3)
    uvec(n_z) = mu_cl - ec0 - kbt * log(ndop_cm3 * 1.0e6_dp / nc_m3)
    do ip = 2, n_z - 1
      uvec(ip) = uvec(1) + (uvec(n_z) - uvec(1)) * real(ip - 1, dp) &
          & / real(n_z - 1, dp)
    end do

    do it = 1, 500
      do ip = 2, n_z - 1
        pp = nv_m3 * exp(min((uvec(ip) - mu_vl) / kbt, 100.0_dp))
        nn = nc_m3 * exp(min((mu_cl - ec0 - uvec(ip)) / kbt, 100.0_dp))
        rhs(ip) = -((uvec(ip + 1) - 2.0_dp * uvec(ip) + uvec(ip - 1)) * off &
            & - cpois * (pp - nn + dopn(ip)))
        bb(ip) = -2.0_dp * off - cpois * (pp + nn) / kbt
      end do
      ! Thomas algorithm for the tridiagonal Newton system
      cp(2) = off / bb(2)
      dq(2) = rhs(2) / bb(2)
      do ip = 3, n_z - 1
        den = bb(ip) - off * cp(ip - 1)
        cp(ip) = off / den
        dq(ip) = (rhs(ip) - off * dq(ip - 1)) / den
      end do
      du(n_z - 1) = dq(n_z - 1)
      do ip = n_z - 2, 2, -1
        du(ip) = dq(ip) - cp(ip) * du(ip + 1)
      end do
      step = 0.0_dp
      do ip = 2, n_z - 1
        du(ip) = sign(min(abs(du(ip)), 0.05_dp), du(ip))
        uvec(ip) = uvec(ip) + du(ip)
        step = max(step, abs(du(ip)))
      end do
      if (step < 1.0e-10_dp) exit
    end do
  end subroutine solve_poisson

  !-----------------------------------------------------------------------------
  ! Transverse-momentum summed local Green's functions of one band at energy
  ! index ie: tilde G<(z), tilde G>(z) in cm^-3 eV^-1 (G< = i tilde G<,
  ! G> = -i tilde G>).
  !   ib = 1: conduction band chain, ib = 2: valence band chain.
  ! The chain has hard-wall ends, a probe broadening gam to the reservoir of the
  ! band at every site (occupation fermi((E-mu)/kT)), the LO-phonon self-energies
  ! pl/pg and, at the defect cell, the defect self-energy s_l/s_g. The diagonals
  ! of G^R, G^<, G^> are obtained from left/right-connected recursions.
  !-----------------------------------------------------------------------------
  subroutine band_energy(ib, ie, mu, tlo, tgo)
    integer, intent(in) :: ib, ie
    real(dp), intent(in) :: mu
    real(dp), intent(out) :: tlo(n_z), tgo(n_z)

    real(dp) :: ee, fo, fh, t2, wk, ez, eperp, w_site, sdl, sdg
    real(dp) :: hh(n_z), tsl(n_z), tsg(n_z), ag(n_z)
    real(dp) :: tslj(n_z), tsgj(n_z), ebe(n_z)
    real(dp) :: laml(n_z), rhol(n_z), lamg(n_z), rhog(n_z)
    complex(dp) :: sg(n_z), glc(n_z), grc(n_z), den
    logical :: in_band
    integer :: j, n

    ee = e_grid(ie)
    fo = fermi((ee - mu) / kbt)
    fh = fermi(-(ee - mu) / kbt)

    sdl = 0.0_dp
    sdg = 0.0_dp
    if (ib == 1) then
      hh(:) = ec0 + u(:) + 2.0_dp * t_c
      ebe(:) = ec0 + u(:)
      t2 = t_c**2
      wk = pk_c * de / dz_cm
      tsl(:) = pl_c(:, ie)
      tsg(:) = pg_c(:, ie)
      if (ee > eps_dl) then
        sdl = aimag(s_l(1, ie))
        sdg = -aimag(s_g(1, ie))
      end if
    else
      hh(:) = ev0 + u(:) - 2.0_dp * t_v
      ebe(:) = ev0 + u(:)
      t2 = t_v**2
      wk = pk_v * de / dz_cm
      tsl(:) = pl_v(:, ie)
      tsg(:) = pg_v(:, ie)
      if (ee < eps_dl) then
        sdl = aimag(s_l(1, ie))
        sdg = -aimag(s_g(1, ie))
      end if
    end if

    tlo(:) = 0.0_dp
    tgo(:) = 0.0_dp

    do j = 1, n_kt
      eperp = (real(j, dp) - 0.5_dp) * de
      if (ib == 1) then
        ez = ee - eperp
      else
        ez = ee + eperp
      end if

      ! Reservoir coupling only where the longitudinal energy lies inside the
      ! local band of the site (finite-bandwidth contact): no broadening
      ! tails in the gap, gap states only from tunnelling (field-induced tails).
      do n = 1, n_z
        if (ib == 1) then
          in_band = (ez >= ebe(n) .and. ez <= ebe(n) + 4.0_dp * t_c)
        else
          in_band = (ez <= ebe(n) .and. ez >= ebe(n) - 4.0_dp * t_v)
        end if
        if (in_band) then
          tslj(n) = tsl(n) + gam * fo
          tsgj(n) = tsg(n) + gam * fh
        else
          tslj(n) = 0.0_dp
          tsgj(n) = 0.0_dp
        end if
        w_site = 0.5_dp * (tslj(n) + tsgj(n)) + 1.0e-7_dp
        sg(n) = cmplx(0.0_dp, -w_site, dp)
      end do
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
        tlo(n) = tlo(n) + wk * ag(n) * (laml(n) + rhol(n) - tslj(n))
        tgo(n) = tgo(n) + wk * ag(n) * (lamg(n) + rhog(n) - tsgj(n))
      end do
    end do
  end subroutine band_energy

  !-----------------------------------------------------------------------------
  ! NEGF self-consistent (SCBA) cycle of the p-i-n device with the defect:
  !   - bands: real-space chains (see band_energy), dressed by the defect
  !     self-energy (libnegf) and by intraband LO-phonon self-energies
  !   - defect: Green's function renormalised by the Dyson-Keldysh equation
  !     (libnegf), occupation and net rate from libnegf
  ! Each pass solves the bands with the current self-energies, then the defect,
  ! then updates the self-energies. The rate is evaluated at the beginning of
  ! each pass. On return: r_conv = recombination rate at the defect cell
  ! [cm^-3 s^-1], nit = passes used.
  !-----------------------------------------------------------------------------
  subroutine negf_scba(muc, muv, r_conv, nit)
    real(dp), intent(in) :: muc, muv
    real(dp), intent(out) :: r_conv
    integer, intent(out) :: nit

    complex(dp) :: gl(1, n_en), gg(1, n_en)
    real(dp) :: nlc(n_z), ngc(n_z), nlv(n_z), ngv(n_z)
    real(dp) :: ee, r_old, r_new
    integer :: ie, it, ip, im

    s_l(:, :) = (0.0_dp, 0.0_dp)
    s_g(:, :) = (0.0_dp, 0.0_dp)
    s_r(:, :) = (0.0_dp, 0.0_dp)
    pl_c(:, :) = 0.0_dp
    pg_c(:, :) = 0.0_dp
    pl_v(:, :) = 0.0_dp
    pg_v(:, :) = 0.0_dp
    r_old = -1.0_dp
    r_new = 0.0_dp
    nit = max_it

    do it = 1, max_it

      ! -- band Green's functions with the current self-energies ------------
      tlc(:, :) = 0.0_dp
      tgc(:, :) = 0.0_dp
      tlv(:, :) = 0.0_dp
      tgv(:, :) = 0.0_dp
      do ie = 1, n_en
        ee = e_grid(ie)
        if (ee >= ec_lo) call band_energy(1, ie, muc, tlc(:, ie), tgc(:, ie))
        if (ee <= ev_hi) call band_energy(2, ie, muv, tlv(:, ie), tgv(:, ie))
        gl(1, ie) = cmplx(0.0_dp, tlc(i_def, ie) + tlv(i_def, ie), dp)
        gg(1, ie) = cmplx(0.0_dp, -(tgc(i_def, ie) + tgv(i_def, ie)), dp)
      end do

      ! -- defect occupation, net rate and defect self-energy (libnegf) -----
      call solve_defect_dyson_keldysh(ens%defects(1), e_grid, w_quad, gl, gg)
      r_new = rho_d * ens%defects(1)%rate_net_el / hbar_evs

      s_l(:, :) = (0.0_dp, 0.0_dp)
      s_g(:, :) = (0.0_dp, 0.0_dp)
      s_r(:, :) = (0.0_dp, 0.0_dp)
      call compute_extended_defect_selfenergy(ens, e_grid, w_quad, s_l, s_g, s_r)

      ! -- LO-phonon self-energies from the local G (mixed) -----------------
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
        pl_c(:, ie) = (1.0_dp - mix) * pl_c(:, ie) + mix * d2_c * nlc(:)
        pg_c(:, ie) = (1.0_dp - mix) * pg_c(:, ie) + mix * d2_c * ngc(:)
        pl_v(:, ie) = (1.0_dp - mix) * pl_v(:, ie) + mix * d2_v * nlv(:)
        pg_v(:, ie) = (1.0_dp - mix) * pg_v(:, ie) + mix * d2_v * ngv(:)
      end do

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

end program pin_defect_current
