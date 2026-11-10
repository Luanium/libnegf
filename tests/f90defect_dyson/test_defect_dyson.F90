!> Unit test for Task 4: solve_defect_dyson_keldysh (Mode B).
!!
!! Same 7-point grid as Task 3:
!!   e = -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4 eV  (de = 0.1 eV)
!!   eps_d = 0.1 eV, e_shift = 0.2 eV
!!   CB resonance: eps_d + e_shift = 0.3 = e_grid(6)  (on-grid, exact hat)
!!   VB resonance: eps_d - e_shift = -0.1 = e_grid(2) (on-grid, exact hat)
!!
!! Physical coupling: m0 = 1.0e-6 so gamma_tot << eta_d = 1.0e-4 eV.
!! This ensures the Mode A / Mode B convergence test is physically meaningful.
program test_defect_dyson

  use ln_defect_srh
  use ln_defect_types
  use ln_precision, only : dp

  implicit none

  integer, parameter :: n_en = 7
  real(dp), parameter :: de = 0.1_dp
  real(dp), parameter :: e0 = -0.2_dp
  real(dp), parameter :: eps_d = 0.1_dp
  real(dp), parameter :: e_shift = 0.2_dp
  real(dp), parameter :: kbt = 0.02585_dp

  real(dp) :: e_grid(n_en), w_quad(n_en)
  integer :: i

  do i = 1, n_en
    e_grid(i) = e0 + real(i - 1, dp) * de
    w_quad(i) = de
  end do

  ! ---------------------------------------------------------------------------
  ! Test 4a: Narrow broadening convergence to Mode A.
  !
  ! Use physical (small) coupling m0 = 1e-6 so gamma_tot << eta_d = 1e-4.
  ! Mode B with a narrow Lorentzian approaches the Mode A delta-function result.
  ! ---------------------------------------------------------------------------
  block
    type(TDefectSpec) :: spec
    complex(dp) :: g_less_diag(1, n_en)
    complex(dp) :: g_great_diag(1, n_en)
    real(dp) :: f_i, B
    real(dp) :: f_occ_A, rate_net_el_A
    real(dp) :: f_occ_B, rate_net_el_B
    real(dp) :: rel_err_focc, rel_err_rate, rel_err_cons
    real(dp) :: f_fd_d, mu_eq

    B = exp(-e_shift / kbt)

    ! Physical weak coupling: gamma_tot ~ 2*m0 << eta_d = 1e-4
    call init_defect_spec(spec, 1, [1], [1.0_dp], 2)
    spec%eps_d = eps_d
    spec%kbt = kbt
    spec%density_weight = 1.0_dp
    spec%ec_cut = eps_d + e_shift - 0.05_dp
    spec%ev_cut = eps_d - e_shift + 0.05_dp
    spec%e_shift(1)   = e_shift       ; spec%e_shift(2)   = e_shift
    spec%m_dc_capt(1) = 1.0e-6_dp     ; spec%m_dc_capt(2) = 0.0_dp
    spec%m_dc_em(1)   = 1.0e-6_dp * B ; spec%m_dc_em(2)   = 0.0_dp
    spec%m_dv_em(1)   = 0.0_dp        ; spec%m_dv_em(2)   = 1.0e-6_dp
    spec%m_dv_capt(1) = 0.0_dp        ; spec%m_dv_capt(2) = 1.0e-6_dp * B
    spec%kernel_type  = DEFECT_KERNEL_CUSTOM

    ! Mode B parameters: narrow window so Mode B → Mode A.
    ! omega_halfwidth must fully contain the Lorentzian (width ~ eta_d + 0.5*gamma)
    ! so we use 100*eta_d to ensure the Lorentzian peak at eps_d is well captured.
    spec%eta_d           = 1.0e-4_dp
    spec%n_omega         = 201
    spec%omega_halfwidth = 100.0_dp * spec%eta_d  ! 0.01 eV — covers Lorentzian

    ! Non-equilibrium Green's functions (same as Task 3b)
    do i = 1, n_en
      if (e_grid(i) >= 0.0_dp) then
        f_i = 1.0_dp / (1.0_dp + exp((e_grid(i) - 0.8_dp) / kbt))
      else
        f_i = 1.0_dp / (1.0_dp + exp((e_grid(i) - (-0.5_dp)) / kbt))
      end if
      g_less_diag(1, i)  = cmplx(0.0_dp,  f_i,            dp)
      g_great_diag(1, i) = cmplx(0.0_dp, -(1.0_dp - f_i), dp)
    end do

    ! Mode A result
    spec%mode = DEFECT_MODE_DELTA
    call solve_defect_occupation(spec, e_grid, w_quad, &
        & g_less_diag, g_great_diag)
    f_occ_A       = spec%f_occ
    rate_net_el_A = spec%rate_net_el

    ! Mode B result
    spec%mode = DEFECT_MODE_DYSON
    call solve_defect_dyson_keldysh(spec, e_grid, w_quad, &
        & g_less_diag, g_great_diag)
    f_occ_B       = spec%f_occ
    rate_net_el_B = spec%rate_net_el

    ! Check 1: f_occ_B ~ f_occ_A within 0.1%
    rel_err_focc = abs(f_occ_B - f_occ_A) / (abs(f_occ_A) + 1.0e-30_dp)
    if (rel_err_focc >= 1.0e-3_dp) then
      write(*, '(a,es14.6,a,es14.6,a,es10.3)') &
          & 'FAIL 4a-1: f_occ_B=', f_occ_B, ' f_occ_A=', f_occ_A, &
          & ' rel_err=', rel_err_focc
      error stop 1
    end if

    ! Check 2: rate_net_el_B ~ rate_net_el_A within 0.1%
    rel_err_rate = abs(rate_net_el_B - rate_net_el_A) &
        & / (abs(rate_net_el_A) + 1.0e-30_dp)
    if (rel_err_rate >= 1.0e-3_dp) then
      write(*, '(a,es14.6,a,es14.6,a,es10.3)') &
          & 'FAIL 4a-2: rate_net_el_B=', rate_net_el_B, &
          & ' rate_net_el_A=', rate_net_el_A, &
          & ' rel_err=', rel_err_rate
      error stop 1
    end if

    ! Check 3: Mode B rate balance
    rel_err_cons = abs(spec%rate_net_el - spec%rate_net_hl) &
        & / (abs(spec%rate_net_el) + 1.0e-30_dp)
    if (rel_err_cons >= 1.0e-10_dp) then
      write(*, '(a,es14.6,a,es14.6,a,es10.3)') &
          & 'FAIL 4a-3: rate_net_el=', spec%rate_net_el, &
          & ' rate_net_hl=', spec%rate_net_hl, &
          & ' rel_err=', rel_err_cons
      error stop 1
    end if

    ! Check 4: thermal equilibrium Mode B => f_occ = f_FD(eps_d)
    mu_eq = 0.05_dp
    do i = 1, n_en
      f_i = 1.0_dp / (1.0_dp + exp((e_grid(i) - mu_eq) / kbt))
      g_less_diag(1, i)  = cmplx(0.0_dp,  f_i,            dp)
      g_great_diag(1, i) = cmplx(0.0_dp, -(1.0_dp - f_i), dp)
    end do
    spec%mode = DEFECT_MODE_DYSON
    call solve_defect_dyson_keldysh(spec, e_grid, w_quad, &
        & g_less_diag, g_great_diag)
    f_fd_d = 1.0_dp / (1.0_dp + exp((eps_d - mu_eq) / kbt))
    if (abs(spec%f_occ - f_fd_d) >= 5.0e-3_dp) then
      write(*, '(a,es14.6,a,es14.6,a,es10.3)') &
          & 'FAIL 4a-4 eq: f_occ=', spec%f_occ, &
          & ' f_FD=', f_fd_d, &
          & ' err=', abs(spec%f_occ - f_fd_d)
      error stop 1
    end if

    call destruct_defect_spec(spec)
  end block

  ! ---------------------------------------------------------------------------
  ! Test 4b: Spectral identities and extended-state conservation in Mode B.
  ! ---------------------------------------------------------------------------
  block
    type(TDefectSpec) :: spec
    type(TDefectEnsemble) :: ensemble
    complex(dp) :: g_less_diag(1, n_en)
    complex(dp) :: g_great_diag(1, n_en)
    complex(dp) :: sigma_less_diag(1, n_en)
    complex(dp) :: sigma_great_diag(1, n_en)
    complex(dp) :: sigma_ret_diag(1, n_en)
    real(dp) :: f_i, B, rel_err_cons
    real(dp) :: pi, sum_rule, keldysh_err, keldysh_max
    real(dp) :: a_d_q, tg_d_less_q, tg_d_great_q, lhs
    real(dp) :: rate_ext_cb, rate_ext_total, tg_less_i, tg_great_i
    integer :: q, ii

    pi = acos(-1.0_dp)
    B = exp(-e_shift / kbt)

    call init_defect_spec(spec, 1, [1], [1.0_dp], 2)
    spec%eps_d = eps_d
    spec%kbt = kbt
    spec%density_weight = 1.0_dp
    spec%ec_cut = eps_d + e_shift - 0.05_dp
    spec%ev_cut = eps_d - e_shift + 0.05_dp
    spec%e_shift(1)   = e_shift       ; spec%e_shift(2)   = e_shift
    spec%m_dc_capt(1) = 1.0e-6_dp     ; spec%m_dc_capt(2) = 0.0_dp
    spec%m_dc_em(1)   = 1.0e-6_dp * B ; spec%m_dc_em(2)   = 0.0_dp
    spec%m_dv_em(1)   = 0.0_dp        ; spec%m_dv_em(2)   = 1.0e-6_dp
    spec%m_dv_capt(1) = 0.0_dp        ; spec%m_dv_capt(2) = 1.0e-6_dp * B
    spec%kernel_type  = DEFECT_KERNEL_CUSTOM
    spec%mode         = DEFECT_MODE_DYSON

    spec%eta_d           = 5.0e-3_dp
    spec%n_omega         = 101
    spec%omega_halfwidth = 10.0_dp * spec%eta_d

    ! Non-equilibrium Green's functions
    do i = 1, n_en
      if (e_grid(i) >= 0.0_dp) then
        f_i = 1.0_dp / (1.0_dp + exp((e_grid(i) - 0.8_dp) / kbt))
      else
        f_i = 1.0_dp / (1.0_dp + exp((e_grid(i) - (-0.5_dp)) / kbt))
      end if
      g_less_diag(1, i)  = cmplx(0.0_dp,  f_i,            dp)
      g_great_diag(1, i) = cmplx(0.0_dp, -(1.0_dp - f_i), dp)
    end do

    call solve_defect_dyson_keldysh(spec, e_grid, w_quad, &
        & g_less_diag, g_great_diag)

    ! Check 1: rate balance
    rel_err_cons = abs(spec%rate_net_el - spec%rate_net_hl) &
        & / (abs(spec%rate_net_el) + 1.0e-30_dp)
    if (rel_err_cons >= 1.0e-10_dp) then
      write(*, '(a,es14.6,a,es14.6,a,es10.3)') &
          & 'FAIL 4b-1: net_el=', spec%rate_net_el, &
          & ' net_hl=', spec%rate_net_hl, ' rel_err=', rel_err_cons
      error stop 1
    end if

    ! Check 2: f_occ in [0,1]
    if (spec%f_occ < 0.0_dp .or. spec%f_occ > 1.0_dp) then
      write(*, '(a,es14.6)') 'FAIL 4b-2: f_occ out of [0,1]:', spec%f_occ
      error stop 1
    end if

    ! Check 3: spectral sum rule int A_d dω/2π = 1
    sum_rule = 0.0_dp
    do q = 1, spec%n_omega
      sum_rule = sum_rule &
          & + (spec%omega_weights(q) / (2.0_dp * pi)) &
          & * (-2.0_dp * aimag(spec%g_d_ret(q)))
    end do
    if (abs(sum_rule - 1.0_dp) >= 1.0e-12_dp) then
      write(*, '(a,es14.6)') 'FAIL 4b-3: spectral sum rule err=', &
          & abs(sum_rule - 1.0_dp)
      error stop 1
    end if

    ! Check 4: Keldysh identity tilde_G_d< + tilde_G_d> = A_d at every q
    keldysh_max = 0.0_dp
    do q = 1, spec%n_omega
      a_d_q        = -2.0_dp * aimag(spec%g_d_ret(q))
      tg_d_less_q  =  aimag(spec%g_d_less(q))
      tg_d_great_q = -aimag(spec%g_d_great(q))
      lhs = tg_d_less_q + tg_d_great_q
      keldysh_err = abs(lhs - a_d_q)
      if (keldysh_err > keldysh_max) keldysh_max = keldysh_err
    end do
    if (keldysh_max >= 1.0e-14_dp) then
      write(*, '(a,es14.6)') 'FAIL 4b-4: Keldysh identity max err=', keldysh_max
      error stop 1
    end if

    ! Check 5: causality Im(G_d_ret(q)) <= 0 for all q
    do q = 1, spec%n_omega
      if (aimag(spec%g_d_ret(q)) > 1.0e-15_dp) then
        write(*, '(a,i0,a,es14.6)') &
            & 'FAIL 4b-5: Im(G_d_ret) > 0 at q=', q, &
            & ' val=', aimag(spec%g_d_ret(q))
        error stop 1
      end if
    end do

    ! Check 6: extended-state particle/current conservation
    call init_defect_ensemble(ensemble, 1, 1)
    ensemble%defects(1) = spec

    sigma_less_diag(:, :)  = (0.0_dp, 0.0_dp)
    sigma_great_diag(:, :) = (0.0_dp, 0.0_dp)
    sigma_ret_diag(:, :)   = (0.0_dp, 0.0_dp)

    call compute_extended_defect_selfenergy(ensemble, e_grid, w_quad, &
        & sigma_less_diag, sigma_great_diag, sigma_ret_diag)

    ! rate_ext_cb = sum over (e_grid(i) >= ec_cut) of collision integral
    ! rate_ext_total = sum over all i
    rate_ext_cb    = 0.0_dp
    rate_ext_total = 0.0_dp
    do ii = 1, n_en
      tg_less_i  =  aimag(g_less_diag(1, ii))
      tg_great_i = -aimag(g_great_diag(1, ii))
      rate_ext_total = rate_ext_total + w_quad(ii) * ( &
          & (-aimag(sigma_great_diag(1, ii))) * tg_less_i &
          & - aimag(sigma_less_diag(1, ii)) * tg_great_i)
      if (e_grid(ii) > spec%eps_d) then
        rate_ext_cb = rate_ext_cb + w_quad(ii) * ( &
            & (-aimag(sigma_great_diag(1, ii))) * tg_less_i &
            & - aimag(sigma_less_diag(1, ii)) * tg_great_i)
      end if
    end do

    if (abs(rate_ext_cb - spec%density_weight * spec%rate_net_el) &
        & / (abs(spec%rate_net_el) + 1.0e-30_dp) >= 1.0e-12_dp) then
      write(*, '(a,es14.6,a,es14.6)') &
          & 'FAIL 4b-6a: rate_ext_cb=', rate_ext_cb, &
          & ' Nd*rate_net_el=', spec%density_weight * spec%rate_net_el
      error stop 1
    end if

    if (abs(rate_ext_total) / (abs(spec%rate_net_el) + 1.0e-30_dp) &
        & >= 1.0e-12_dp) then
      write(*, '(a,es14.6)') &
          & 'FAIL 4b-6b: rate_ext_total should be zero, got', rate_ext_total
      error stop 1
    end if

    call destruct_defect_ensemble(ensemble)
    call destruct_defect_spec(spec)
  end block

  write(*, '(a)') 'Task 4: test_defect_dyson PASSED'

end program test_defect_dyson
