!> Unit test for Task 3: ln_defect_srh routines.
!! Tests solve_defect_occupation and compute_extended_defect_selfenergy.
!!
!! Strategy: 7-point uniform grid with resonance energies placed EXACTLY
!! on grid points so hat-interpolation is exact to machine precision.
!!
!! Grid: e = -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4 eV  (de = 0.1 eV)
!! eps_d = 0.1 eV, e_shift = 0.2 eV
!! CB resonance: eps_d + e_shift = 0.3 = e_grid(6)
!! VB resonance: eps_d - e_shift = -0.1 = e_grid(2)
program test_defect_srh_unit

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
  real(dp), parameter :: mu_eq = 0.05_dp

  real(dp) :: e_grid(n_en), w_quad(n_en)
  integer :: i

  do i = 1, n_en
    e_grid(i) = e0 + real(i - 1, dp) * de
    w_quad(i) = de
  end do
  ! e_grid(2) = -0.1 (VB resonance), e_grid(6) = 0.3 (CB resonance)

  ! Helper: build a standard 2-mode custom-kernel spec
  ! (mode 1 = CB resonance at eps_d+e_shift, mode 2 = VB resonance at eps_d-e_shift)
  ! ec_cut set just below CB resonance, ev_cut just above VB resonance.

  ! ---------------------------------------------------------------------------
  ! Test 3a: Thermal equilibrium => f_occ = f_FD(eps_d), zero net rates.
  ! ---------------------------------------------------------------------------
  block
    type(TDefectSpec) :: spec
    complex(dp) :: g_less_diag(1, n_en)
    complex(dp) :: g_great_diag(1, n_en)
    real(dp) :: f_i, f_fd_d, B, err_focc, err_el, err_hl

    B = exp(-e_shift / kbt)   ! ~4.35e-4, well above underflow

    call init_defect_spec(spec, 1, [1], [1.0_dp], 2)
    spec%eps_d = eps_d
    spec%kbt = kbt
    spec%density_weight = 1.0_dp
    spec%ec_cut = eps_d + e_shift - 0.05_dp   ! = 0.25: CB resonance 0.3 passes
    spec%ev_cut = eps_d - e_shift + 0.05_dp   ! = -0.05: VB resonance -0.1 passes
    ! Mode 1: CB — m_dc_capt=1, m_dc_em=B, m_dv_*=0
    ! Mode 2: VB — m_dv_em=1, m_dv_capt=B, m_dc_*=0
    spec%e_shift(1)   = e_shift  ; spec%e_shift(2)   = e_shift
    spec%m_dc_capt(1) = 1.0_dp   ; spec%m_dc_capt(2) = 0.0_dp
    spec%m_dc_em(1)   = B        ; spec%m_dc_em(2)   = 0.0_dp
    spec%m_dv_em(1)   = 0.0_dp   ; spec%m_dv_em(2)   = 1.0_dp
    spec%m_dv_capt(1) = 0.0_dp   ; spec%m_dv_capt(2) = B
    spec%kernel_type  = DEFECT_KERNEL_CUSTOM

    ! Equilibrium Green's functions
    do i = 1, n_en
      f_i = 1.0_dp / (1.0_dp + exp((e_grid(i) - mu_eq) / kbt))
      g_less_diag(1, i)  = cmplx(0.0_dp,  f_i,            dp)
      g_great_diag(1, i) = cmplx(0.0_dp, -(1.0_dp - f_i), dp)
    end do

    call solve_defect_occupation(spec, e_grid, w_quad, &
        & g_less_diag, g_great_diag)

    ! Equilibrium identity: f_occ = f_FD(eps_d)
    ! Proof: with on-grid hat, tilde_Sigma_dc< = M_capt*f_c, tilde_Sigma_dc> = M_em*(1-f_c) = B*(1-f_c).
    ! Similarly tilde_Sigma_dv< = B*f_v, tilde_Sigma_dv> = (1-f_v).
    ! f_occ = (f_c + B*f_v) / (f_c + B*(1-f_c) + B*f_v + (1-f_v))
    ! At equilibrium with detailed balance B = exp(-e_shift/kbt):
    !   f_c = f_FD(e_c) = f_FD(eps_d + e_shift)
    !   It follows algebraically that f_occ = f_FD(eps_d).
    f_fd_d = 1.0_dp / (1.0_dp + exp((eps_d - mu_eq) / kbt))
    err_focc = abs(spec%f_occ - f_fd_d)
    if (err_focc >= 1.0e-12_dp) then
      write(*, '(a,es16.8,a,es16.8,a,es10.3)') &
          & 'FAIL 3a f_occ: got', spec%f_occ, ' expected', f_fd_d, ' err', err_focc
      error stop 1
    end if

    ! Zero net rates at equilibrium
    err_el = abs(spec%rate_net_el)
    err_hl = abs(spec%rate_net_hl)
    if (err_el >= 1.0e-13_dp .or. err_hl >= 1.0e-13_dp) then
      write(*, '(a,es14.6,a,es14.6)') &
          & 'FAIL 3a rates: net_el=', spec%rate_net_el, ' net_hl=', spec%rate_net_hl
      error stop 1
    end if

    call destruct_defect_spec(spec)
  end block

  ! ---------------------------------------------------------------------------
  ! Test 3b: Non-equilibrium => rate balance, causal sigma_ret, and exact
  ! extended-state particle conservation (Verify 3).
  ! ---------------------------------------------------------------------------
  block
    type(TDefectSpec) :: spec
    type(TDefectEnsemble) :: ensemble
    complex(dp) :: g_less_diag(1, n_en)
    complex(dp) :: g_great_diag(1, n_en)
    complex(dp) :: sigma_less_diag(1, n_en)
    complex(dp) :: sigma_great_diag(1, n_en)
    complex(dp) :: sigma_ret_diag(1, n_en)
    real(dp) :: f_i, B, r_el, r_hl, rel_err
    real(dp) :: sret_expect, fdt_err, fdt_max
    real(dp) :: rate_ext_cb, rate_ext_total, tg_less_i, tg_great_i
    integer :: ii

    B = exp(-e_shift / kbt)

    call init_defect_spec(spec, 1, [1], [1.0_dp], 2)
    spec%eps_d = eps_d
    spec%kbt = kbt
    spec%density_weight = 1.0_dp
    spec%ec_cut = eps_d + e_shift - 0.05_dp
    spec%ev_cut = eps_d - e_shift + 0.05_dp
    spec%e_shift(1)   = e_shift  ; spec%e_shift(2)   = e_shift
    spec%m_dc_capt(1) = 1.0_dp   ; spec%m_dc_capt(2) = 0.0_dp
    spec%m_dc_em(1)   = B        ; spec%m_dc_em(2)   = 0.0_dp
    spec%m_dv_em(1)   = 0.0_dp   ; spec%m_dv_em(2)   = 1.0_dp
    spec%m_dv_capt(1) = 0.0_dp   ; spec%m_dv_capt(2) = B
    spec%kernel_type  = DEFECT_KERNEL_CUSTOM

    ! Non-equilibrium: upper half uses mu_c, lower half uses mu_v
    do i = 1, n_en
      if (e_grid(i) >= 0.0_dp) then
        f_i = 1.0_dp / (1.0_dp + exp((e_grid(i) - 0.8_dp) / kbt))
      else
        f_i = 1.0_dp / (1.0_dp + exp((e_grid(i) - (-0.5_dp)) / kbt))
      end if
      g_less_diag(1, i)  = cmplx(0.0_dp,  f_i,            dp)
      g_great_diag(1, i) = cmplx(0.0_dp, -(1.0_dp - f_i), dp)
    end do

    call solve_defect_occupation(spec, e_grid, w_quad, &
        & g_less_diag, g_great_diag)

    ! Verify 1: rate_net_el == rate_net_hl (algebraic identity)
    r_el = spec%rate_net_el
    r_hl = spec%rate_net_hl
    rel_err = abs(r_el - r_hl) / (abs(r_el) + 1.0e-30_dp)
    if (rel_err >= 1.0e-12_dp) then
      write(*, '(a,es14.6,a,es14.6,a,es10.3)') &
          & 'FAIL 3b-1: net_el=', r_el, ' net_hl=', r_hl, ' rel_err=', rel_err
      error stop 1
    end if

    ! Verify 1b: non-zero recombination (out-of-equilibrium system)
    if (abs(r_el) < 1.0e-20_dp) then
      write(*, '(a,es14.6)') 'FAIL 3b-1b: rate_net_el should be nonzero, got', r_el
      error stop 1
    end if

    ! Compute extended self-energies
    call init_defect_ensemble(ensemble, 1, 1)
    ensemble%defects(1) = spec

    sigma_less_diag(:, :)  = (0.0_dp, 0.0_dp)
    sigma_great_diag(:, :) = (0.0_dp, 0.0_dp)
    sigma_ret_diag(:, :)   = (0.0_dp, 0.0_dp)

    call compute_extended_defect_selfenergy(ensemble, e_grid, w_quad, &
        & sigma_less_diag, sigma_great_diag, sigma_ret_diag)

    ! Verify 2: causal retarded self-energy
    ! Im(Sigma^R) = +0.5*(Im(Sigma>) - Im(Sigma<)) <= 0
    fdt_max = 0.0_dp
    do ii = 1, n_en
      sret_expect = 0.5_dp * &
          & (aimag(sigma_great_diag(1, ii)) - aimag(sigma_less_diag(1, ii)))
      fdt_err = abs(aimag(sigma_ret_diag(1, ii)) - sret_expect)
      if (fdt_err > fdt_max) fdt_max = fdt_err
      ! Causality: Im(Sigma^R) <= 0
      if (aimag(sigma_ret_diag(1, ii)) > 1.0e-15_dp) then
        write(*, '(a,i0,a,es14.6)') &
            & 'FAIL 3b-2 causality at i=', ii, &
            & ': Im(Sigma_ret)=', aimag(sigma_ret_diag(1, ii))
        error stop 1
      end if
    end do
    if (fdt_max >= 1.0e-14_dp) then
      write(*, '(a,es14.6)') 'FAIL 3b-2 FDT max error=', fdt_max
      error stop 1
    end if

    ! Verify 3: exact extended-state particle conservation.
    ! The collision integral sum_i w_quad(i) * (Sigma>(i)*G<(i) - Sigma<(i)*G>(i))
    ! equals spec%density_weight * spec%rate_net_el when restricted to CB energies,
    ! and the TOTAL integral over all energies equals zero (equal in/out rates).
    !
    ! In tilde notation:
    !   integrand(i) = tilde_Sigma>(i) * tilde_G<(i) - tilde_Sigma<(i) * tilde_G>(i)
    !                = (-Im(Sigma>(i))) * Im(G<(i)) - Im(Sigma<(i)) * (-Im(G>(i)))
    rate_ext_cb    = 0.0_dp
    rate_ext_total = 0.0_dp
    do ii = 1, n_en
      tg_less_i  =  aimag(g_less_diag(1, ii))    ! Im(G<) = tilde_G<
      tg_great_i = -aimag(g_great_diag(1, ii))   ! -Im(G>) = tilde_G>
      rate_ext_total = rate_ext_total + w_quad(ii) * ( &
          & (-aimag(sigma_great_diag(1, ii))) * tg_less_i &
          & - aimag(sigma_less_diag(1, ii)) * tg_great_i)
      ! CB band: e_grid(i) >= ec_cut
      if (e_grid(ii) >= spec%ec_cut) then
        rate_ext_cb = rate_ext_cb + w_quad(ii) * ( &
            & (-aimag(sigma_great_diag(1, ii))) * tg_less_i &
            & - aimag(sigma_less_diag(1, ii)) * tg_great_i)
      end if
    end do

    if (abs(rate_ext_cb - spec%density_weight * spec%rate_net_el) >= 1.0e-12_dp) then
      write(*, '(a,es14.6,a,es14.6)') &
          & 'FAIL 3b-3a: rate_ext_cb=', rate_ext_cb, &
          & ' Nd*rate_net_el=', spec%density_weight * spec%rate_net_el
      error stop 1
    end if
    if (abs(rate_ext_total) >= 1.0e-12_dp) then
      write(*, '(a,es14.6)') 'FAIL 3b-3b: rate_ext_total should be zero, got', rate_ext_total
      error stop 1
    end if

    call destruct_defect_ensemble(ensemble)
    call destruct_defect_spec(spec)
  end block

  ! ---------------------------------------------------------------------------
  ! Test 3c: ec_cut guard in compute_extended_defect_selfenergy.
  ! Set ec_cut ABOVE the CB resonance so the CB mode is suppressed.
  ! Then rate_ext_total must still be zero (no mismatch from missing guard).
  ! ---------------------------------------------------------------------------
  block
    type(TDefectSpec) :: spec
    type(TDefectEnsemble) :: ensemble
    complex(dp) :: g_less_diag(1, n_en)
    complex(dp) :: g_great_diag(1, n_en)
    complex(dp) :: sigma_less_diag(1, n_en)
    complex(dp) :: sigma_great_diag(1, n_en)
    complex(dp) :: sigma_ret_diag(1, n_en)
    real(dp) :: f_i, B, rate_ext_total, tg_less_i, tg_great_i
    integer :: ii

    B = exp(-e_shift / kbt)

    call init_defect_spec(spec, 1, [1], [1.0_dp], 2)
    spec%eps_d = eps_d
    spec%kbt = kbt
    spec%density_weight = 1.0_dp
    ! ec_cut ABOVE CB resonance (0.3): CB mode suppressed in both routines
    spec%ec_cut = eps_d + e_shift + 0.05_dp   ! = 0.35 > 0.3
    spec%ev_cut = eps_d - e_shift + 0.05_dp   ! = -0.05: VB mode active
    spec%e_shift(1)   = e_shift  ; spec%e_shift(2)   = e_shift
    spec%m_dc_capt(1) = 1.0_dp   ; spec%m_dc_capt(2) = 0.0_dp
    spec%m_dc_em(1)   = B        ; spec%m_dc_em(2)   = 0.0_dp
    spec%m_dv_em(1)   = 0.0_dp   ; spec%m_dv_em(2)   = 1.0_dp
    spec%m_dv_capt(1) = 0.0_dp   ; spec%m_dv_capt(2) = B
    spec%kernel_type  = DEFECT_KERNEL_CUSTOM

    do i = 1, n_en
      if (e_grid(i) >= 0.0_dp) then
        f_i = 1.0_dp / (1.0_dp + exp((e_grid(i) - 0.8_dp) / kbt))
      else
        f_i = 1.0_dp / (1.0_dp + exp((e_grid(i) - (-0.5_dp)) / kbt))
      end if
      g_less_diag(1, i)  = cmplx(0.0_dp,  f_i,            dp)
      g_great_diag(1, i) = cmplx(0.0_dp, -(1.0_dp - f_i), dp)
    end do

    call solve_defect_occupation(spec, e_grid, w_quad, &
        & g_less_diag, g_great_diag)

    call init_defect_ensemble(ensemble, 1, 1)
    ensemble%defects(1) = spec

    sigma_less_diag(:, :)  = (0.0_dp, 0.0_dp)
    sigma_great_diag(:, :) = (0.0_dp, 0.0_dp)
    sigma_ret_diag(:, :)   = (0.0_dp, 0.0_dp)

    call compute_extended_defect_selfenergy(ensemble, e_grid, w_quad, &
        & sigma_less_diag, sigma_great_diag, sigma_ret_diag)

    ! Only VB mode contributes. Total collision integral must still be zero
    ! because sigma only has non-zero entries at VB resonance and f_occ was
    ! solved consistently with the same guard.
    rate_ext_total = 0.0_dp
    do ii = 1, n_en
      tg_less_i  =  aimag(g_less_diag(1, ii))
      tg_great_i = -aimag(g_great_diag(1, ii))
      rate_ext_total = rate_ext_total + w_quad(ii) * ( &
          & (-aimag(sigma_great_diag(1, ii))) * tg_less_i &
          & - aimag(sigma_less_diag(1, ii)) * tg_great_i)
    end do
    if (abs(rate_ext_total) >= 1.0e-12_dp) then
      write(*, '(a,es14.6)') &
          & 'FAIL 3c: rate_ext_total with suppressed CB mode =', rate_ext_total
      error stop 1
    end if

    call destruct_defect_ensemble(ensemble)
    call destruct_defect_spec(spec)
  end block

  write(*, '(a)') 'Task 3: test_defect_srh_unit PASSED'

end program test_defect_srh_unit
