!> Unit test for Task 5: SRH defect SCBA update.
!!
!! Tests the fixed-point SCBA iteration for the SRH defect self-energy on a
!! minimal two-band model:
!!   - 8-point uniform energy grid covering [-0.4, 0.3] eV (de = 0.1 eV)
!!   - Two orbitals: orbital 1 = VB, orbital 2 = CB
!!   - Midgap defect at eps_d = 0.0 eV, e_shift = 0.2 eV
!!   - CB resonance at eps_d + e_shift = 0.2 eV (inside grid)
!!   - VB resonance at eps_d - e_shift = -0.2 eV (inside grid)
!!   - Non-equilibrium G</G>: VB filled (mu_v = -0.3 eV), CB empty (mu_c = 0.5 eV)
program test_defect_scba

  use ln_defect_srh
  use ln_defect_types
  use ln_defect_kernel
  use ln_precision, only : dp
  use libnegf, only : TNegf, init_negf, destroy_negf, &
      & set_defect_srh_mpr, get_defect_srh_rates

  implicit none

  integer, parameter :: n_en   = 8
  integer, parameter :: n_orb  = 2
  real(dp), parameter :: de    = 0.1_dp
  real(dp), parameter :: e0    = -0.4_dp
  real(dp), parameter :: eps_d = 0.0_dp
  real(dp), parameter :: e_shift    = 0.2_dp
  real(dp), parameter :: kbt        = 0.02585_dp
  real(dp), parameter :: m0_dc      = 1.0e-5_dp
  real(dp), parameter :: m0_dv      = 1.0e-5_dp
  real(dp), parameter :: hbar_omega0 = 0.1_dp
  real(dp), parameter :: s_hr       = 1.0_dp
  integer, parameter  :: n_modes    = 5
  integer, parameter  :: max_scba   = 100
  real(dp), parameter :: scba_tol   = 1.0e-9_dp

  real(dp) :: e_grid(n_en), w_quad(n_en)
  integer :: i

  do i = 1, n_en
    e_grid(i) = e0 + real(i - 1, dp) * de
    w_quad(i) = de
  end do

  ! =========================================================================
  ! Step 1: Zero-defect baseline.
  ! With m0=0, sigma</>= 0 and rate_net_el should be exactly 0.
  ! =========================================================================
  block
    type(TDefectSpec) :: spec

    call init_defect_spec(spec, 1, [2], [1.0_dp], n_modes)
    spec%eps_d          = eps_d
    spec%kbt            = kbt
    spec%density_weight = 1.0_dp
    spec%hbar_omega0    = hbar_omega0
    spec%huang_rhys_n   = s_hr
    spec%huang_rhys_p   = s_hr
    spec%m0_dc          = 0.0_dp
    spec%m0_dv          = 0.0_dp
    spec%ec_cut         = eps_d + e_shift - 0.05_dp
    spec%ev_cut         = eps_d - e_shift + 0.05_dp
    spec%mode           = DEFECT_MODE_DELTA
    call build_mpr_kernel(spec)

    block
      complex(dp) :: g_less(1, n_en), g_great(1, n_en)
      real(dp) :: f_i

      do i = 1, n_en
        f_i = 1.0_dp / (1.0_dp + exp((e_grid(i) - 0.4_dp) / kbt))
        g_less(1, i)  = cmplx(0.0_dp,  f_i,            dp)
        g_great(1, i) = cmplx(0.0_dp, -(1.0_dp - f_i), dp)
      end do

      call solve_defect_occupation(spec, e_grid, w_quad, g_less, g_great)

      if (abs(spec%rate_net_el) > 1.0e-30_dp) then
        write(*, '(a,es14.6)') &
            & 'FAIL step1: rate_net_el should be 0 for zero coupling, got', &
            & spec%rate_net_el
        error stop 1
      end if
    end block

    call destruct_defect_spec(spec)
  end block

  ! =========================================================================
  ! Step 2: Midgap MPR defect with two orbitals (VB=1, CB=2).
  ! Run a coupled Dyson-Keldysh SCBA that self-consistently updates both
  ! the defect occupation and the dressed Green functions.
  ! =========================================================================
  block
    type(TDefectEnsemble) :: ensemble

    ! Dressed Green's functions (updated each SCBA iteration)
    complex(dp) :: g_less(n_orb, n_en)
    complex(dp) :: g_great(n_orb, n_en)
    complex(dp) :: g_ret(n_orb, n_en)

    ! Bare (contact-only) quantities
    complex(dp) :: g_ret_0(n_orb, n_en)
    complex(dp) :: sigma_0_less(n_orb, n_en)
    complex(dp) :: sigma_0_great(n_orb, n_en)

    ! Defect self-energies (zeroed before each SCBA iteration)
    complex(dp) :: sigma_less(n_orb, n_en)
    complex(dp) :: sigma_great(n_orb, n_en)
    complex(dp) :: sigma_ret(n_orb, n_en)
    complex(dp) :: sigma_less_prev(n_orb, n_en)

    real(dp) :: f_vb, f_cb, gamma_0
    real(dp) :: mu_v, mu_c
    real(dp) :: f_occ_prev, scba_err
    real(dp) :: rate_ext_total, rate_ext_cb
    real(dp) :: tg_less_m, tg_great_m
    real(dp) :: abs_g_ret_sq
    integer :: iter, m_orb
    logical :: converged

    ! Contact broadening parameters
    gamma_0 = 0.05_dp
    mu_v    = -0.3_dp
    mu_c    =  0.5_dp

    call init_defect_ensemble(ensemble, 1, 1)

    associate(spec => ensemble%defects(1))
      call init_defect_spec(spec, 2, [1, 2], &
          & [0.5_dp, 0.5_dp], n_modes)
      spec%eps_d          = eps_d
      spec%kbt            = kbt
      spec%density_weight = 1.0_dp
      spec%hbar_omega0    = hbar_omega0
      spec%huang_rhys_n   = s_hr
      spec%huang_rhys_p   = s_hr
      spec%m0_dc          = m0_dc
      spec%m0_dv          = m0_dv
      spec%ec_cut         = eps_d + e_shift - 0.05_dp
      spec%ev_cut         = eps_d - e_shift + 0.05_dp
      spec%mode           = DEFECT_MODE_DELTA
      spec%pl_index       = 1
      call build_mpr_kernel(spec)
    end associate

    ! --- Define bare (contact-only) Green's functions ---
    ! g_ret_0 = 1 / (i * 0.5 * gamma_0) for a featureless broadened level
    do i = 1, n_en
      f_vb = 1.0_dp / (1.0_dp + exp((e_grid(i) - mu_v) / kbt))
      f_cb = 1.0_dp / (1.0_dp + exp((e_grid(i) - mu_c) / kbt))

      ! Bare retarded GF: 1 / (i * gamma_0 / 2) for both orbitals
      g_ret_0(1, i) = 1.0_dp / cmplx(0.0_dp, 0.5_dp * gamma_0, dp)
      g_ret_0(2, i) = 1.0_dp / cmplx(0.0_dp, 0.5_dp * gamma_0, dp)

      ! Contact self-energies (imaginary only)
      sigma_0_less(1, i)  = cmplx(0.0_dp,  gamma_0 * f_vb, dp)
      sigma_0_great(1, i) = cmplx(0.0_dp, -gamma_0 * (1.0_dp - f_vb), dp)
      sigma_0_less(2, i)  = cmplx(0.0_dp,  gamma_0 * f_cb, dp)
      sigma_0_great(2, i) = cmplx(0.0_dp, -gamma_0 * (1.0_dp - f_cb), dp)
    end do

    ! --- Initialise defect self-energies to zero ---
    sigma_less(:, :)  = (0.0_dp, 0.0_dp)
    sigma_great(:, :) = (0.0_dp, 0.0_dp)
    sigma_ret(:, :)   = (0.0_dp, 0.0_dp)

    ! --- Coupled Dyson-Keldysh SCBA loop ---
    f_occ_prev = 0.5_dp
    converged  = .false.
    scba_err   = 1.0_dp
    iter       = 0

    do while (.not. converged .and. iter < max_scba)

      ! Step a: compute dressed GF from bare + defect self-energy
      do i = 1, n_en
        do m_orb = 1, n_orb
          g_ret(m_orb, i) = 1.0_dp / &
              & (1.0_dp / g_ret_0(m_orb, i) - sigma_ret(m_orb, i))
          abs_g_ret_sq = real(g_ret(m_orb, i))**2 + aimag(g_ret(m_orb, i))**2
          g_less(m_orb, i) = abs_g_ret_sq * &
              & (sigma_0_less(m_orb, i) + sigma_less(m_orb, i))
          g_great(m_orb, i) = abs_g_ret_sq * &
              & (sigma_0_great(m_orb, i) + sigma_great(m_orb, i))
        end do
      end do

      ! Step b: solve defect occupation from dressed G</G>
      call solve_defect_occupation(ensemble%defects(1), e_grid, w_quad, &
          & g_less, g_great)

      ! Step c: save previous values for convergence check (before zeroing sigma)
      sigma_less_prev(:, :) = sigma_less(:, :)

      ! Step d: zero out self-energies before recomputing
      sigma_less(:, :)  = (0.0_dp, 0.0_dp)
      sigma_great(:, :) = (0.0_dp, 0.0_dp)
      sigma_ret(:, :)   = (0.0_dp, 0.0_dp)

      ! Step e: accumulate ensemble rates and compute new self-energies
      ensemble%pl_recomb_rate(1) = &
          & ensemble%defects(1)%density_weight * ensemble%defects(1)%rate_net_el
      ensemble%total_recomb_rate = ensemble%pl_recomb_rate(1)

      call compute_extended_defect_selfenergy(ensemble, e_grid, w_quad, &
          & sigma_less, sigma_great, sigma_ret)

      ! Step f: check convergence on f_occ and sigma_less
      ! f_occ_prev was set at the top of the previous iteration; compare with
      ! the new f_occ computed in step b and the new sigma from step e.
      scba_err = abs(ensemble%defects(1)%f_occ - f_occ_prev)
      if (scba_err < scba_tol .and. &
          & maxval(abs(aimag(sigma_less) - aimag(sigma_less_prev))) < scba_tol) then
        converged = .true.
      end if
      f_occ_prev = ensemble%defects(1)%f_occ

      iter = iter + 1
    end do

    ! Step 4: final dressed GF update with converged sigma
    do i = 1, n_en
      do m_orb = 1, n_orb
        g_ret(m_orb, i) = 1.0_dp / &
            & (1.0_dp / g_ret_0(m_orb, i) - sigma_ret(m_orb, i))
        abs_g_ret_sq = real(g_ret(m_orb, i))**2 + aimag(g_ret(m_orb, i))**2
        g_less(m_orb, i) = abs_g_ret_sq * &
            & (sigma_0_less(m_orb, i) + sigma_less(m_orb, i))
        g_great(m_orb, i) = abs_g_ret_sq * &
            & (sigma_0_great(m_orb, i) + sigma_great(m_orb, i))
      end do
    end do

    ! Verify 1: SCBA converged within < 30 iterations
    if (.not. converged) then
      write(*, '(a,i0,a,es10.3)') &
          & 'FAIL step2-1: SCBA did not converge in ', max_scba, &
          & ' iterations; last scba_err=', scba_err
      error stop 1
    end if
    if (iter >= 30) then
      write(*, '(a,i0,a)') &
          & 'FAIL step2-1: SCBA converged in ', iter, &
          & ' iterations (expected < 30)'
      error stop 1
    end if

    ! Verify 2: Finite interband recombination rate > 0
    if (ensemble%defects(1)%rate_net_el <= 0.0_dp) then
      write(*, '(a,es14.6)') &
          & 'FAIL step2-2: rate_net_el should be > 0, got', &
          & ensemble%defects(1)%rate_net_el
      error stop 1
    end if

    ! Verify 3: Current conservation via extended self-energy.
    ! Reuse the converged sigma_less/great (already filled from last iteration).
    !
    ! rate_ext_total = sum_i w_i * sum_m [(-Im Sigma>) Im(G<) - Im(Sigma<) (-Im G>)]
    ! Should be zero (particle conservation).
    ! rate_ext_cb sums ALL orbitals for e_grid(i) > eps_d (CB region).
    rate_ext_total = 0.0_dp
    rate_ext_cb    = 0.0_dp
    do i = 1, n_en
      do m_orb = 1, n_orb
        tg_less_m  =  aimag(g_less(m_orb, i))
        tg_great_m = -aimag(g_great(m_orb, i))
        rate_ext_total = rate_ext_total + w_quad(i) * ( &
            & (-aimag(sigma_great(m_orb, i))) * tg_less_m &
            & - aimag(sigma_less(m_orb, i)) * tg_great_m)
      end do
      ! CB collision integral: sum all orbitals for energies above defect level
      if (e_grid(i) > eps_d) then
        do m_orb = 1, n_orb
          tg_less_m  =  aimag(g_less(m_orb, i))
          tg_great_m = -aimag(g_great(m_orb, i))
          rate_ext_cb = rate_ext_cb + w_quad(i) * ( &
              & (-aimag(sigma_great(m_orb, i))) * tg_less_m &
              & - aimag(sigma_less(m_orb, i)) * tg_great_m)
        end do
      end if
    end do

    ! Total collision integral should be zero (conservation)
    if (abs(rate_ext_total) / &
        & (abs(ensemble%defects(1)%rate_net_el) + 1.0e-30_dp) >= 1.0e-10_dp) then
      write(*, '(a,es14.6,a,es14.6)') &
          & 'FAIL step2-3a: rate_ext_total=', rate_ext_total, &
          & ' (expected 0), rate_net_el=', ensemble%defects(1)%rate_net_el
      error stop 1
    end if

    ! CB collision integral should equal density_weight * rate_net_el
    if (abs(rate_ext_cb &
        & - ensemble%defects(1)%density_weight * ensemble%defects(1)%rate_net_el) &
        & / (abs(ensemble%defects(1)%rate_net_el) + 1.0e-30_dp) >= 1.0e-12_dp) then
      write(*, '(a,es14.6,a,es14.6)') &
          & 'FAIL step2-3b: rate_ext_cb=', rate_ext_cb, &
          & ' Nd*rate_net_el=', &
          & ensemble%defects(1)%density_weight * ensemble%defects(1)%rate_net_el
      error stop 1
    end if

    call destruct_defect_ensemble(ensemble)
  end block

  ! =========================================================================
  ! Step 3: Smoke-test the libnegf public API for SRH defect setup and
  ! rate readout, using the same physical parameters as Step 2.
  ! =========================================================================
  block
    type(TNegf) :: negf_obj

    real(dp) :: f_occ_arr(1), rate_el_arr(1), rate_hl_arr(1)
    real(dp) :: pl_rate(1), total_rate
    integer :: orb_start(1), orb_end(1), pl_idx(1)
    real(dp) :: eps_d_arr(1), dw_arr(1), ec_arr(1), ev_arr(1)

    ! Minimal init: no contacts, 1 principal layer set directly on str
    call init_negf(negf_obj)
    negf_obj%str%num_PLs = 1

    orb_start(1) = 1
    orb_end(1)   = 2
    pl_idx(1)    = 1
    eps_d_arr(1) = eps_d
    dw_arr(1)    = 1.0_dp
    ec_arr(1)    = eps_d + e_shift - 0.05_dp
    ev_arr(1)    = eps_d - e_shift + 0.05_dp

    call set_defect_srh_mpr(negf_obj, 1, pl_idx, orb_start, orb_end, &
        & eps_d_arr, dw_arr, kbt, hbar_omega0, s_hr, s_hr, &
        & m0_dc, m0_dv, n_modes, ec_arr, ev_arr, DEFECT_MODE_DELTA)

    ! Check is_active flag set by set_defect_srh_mpr
    if (.not. (negf_obj%defect_srh%is_active .eqv. .true.)) then
      write(*, '(a)') 'FAIL step3-1: defect_srh%is_active should be .true.'
      error stop 1
    end if

    ! Exercise alloc/dealloc_ensemble_sigma (allocates the sigma cache arrays)
    call alloc_ensemble_sigma(negf_obj%defect_srh, n_orb, n_en)
    if (.not. allocated(negf_obj%defect_srh%sigma_less)) then
      write(*, '(a)') 'FAIL step3-2: sigma_less not allocated after alloc_ensemble_sigma'
      error stop 1
    end if
    call dealloc_ensemble_sigma(negf_obj%defect_srh)
    if (allocated(negf_obj%defect_srh%sigma_less)) then
      write(*, '(a)') 'FAIL step3-3: sigma_less still allocated after dealloc_ensemble_sigma'
      error stop 1
    end if

    ! get_defect_srh_rates: arrays must match num_defects=1, pl_recomb_rate size=1
    call get_defect_srh_rates(negf_obj, f_occ_arr, rate_el_arr, rate_hl_arr, &
        & pl_rate, total_rate)
    if (size(f_occ_arr) /= negf_obj%defect_srh%num_defects) then
      write(*, '(a)') 'FAIL step3-4: f_occ_arr size mismatch'
      error stop 1
    end if

    call destroy_negf(negf_obj)

    write(*, '(a)') 'Step 3 (libnegf API): PASSED'
  end block

  write(*, '(a)') 'Task 5: test_defect_scba PASSED'

end program test_defect_scba
