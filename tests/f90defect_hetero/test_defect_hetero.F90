!> Unit test for Task 6: band-resolved inter-PL currents from SRH rates.
!!
!! Test 6a: digital current-conservation step at a defective heterointerface.
!!   1 defect at pl_index=3 in a 6-PL structure.  After SCBA convergence,
!!   get_band_resolved_pl_currents must reproduce the Aeberhard step profile
!!   (continuity mode) AND must give identical results in spectral mode when
!!   fed the collision integrand from compute_extended_defect_selfenergy.
!!
!! Test 6b: field-enhanced subgap tunneling recombination (Franz-Keldysh,
!!   Aeberhard Fig. 3b).  Increasing the subgap spectral tail (gamma_tail)
!!   to mimic a finite electric field raises the recombination rate.
program test_defect_hetero

  use ln_defect_srh
  use ln_defect_types
  use ln_defect_kernel
  use ln_precision, only : dp

  implicit none

  ! -----------------------------------------------------------------------
  ! Shared energy-grid parameters (identical to Task 5)
  ! -----------------------------------------------------------------------
  integer, parameter :: n_en  = 8
  integer, parameter :: n_orb = 2
  real(dp), parameter :: de     = 0.1_dp
  real(dp), parameter :: e0     = -0.4_dp
  real(dp), parameter :: eps_d  = 0.0_dp
  real(dp), parameter :: e_shift = 0.2_dp
  real(dp), parameter :: kbt    = 0.02585_dp
  real(dp), parameter :: m_dc_capt_val = 1.0e-5_dp
  real(dp), parameter :: m_dv_em_val   = 1.0e-5_dp
  real(dp), parameter :: ec_cut_val    =  0.15_dp
  real(dp), parameter :: ev_cut_val    = -0.15_dp
  real(dp), parameter :: mu_v = -0.3_dp
  real(dp), parameter :: mu_c =  0.5_dp
  real(dp), parameter :: gamma_0 = 0.05_dp
  integer,  parameter :: max_scba = 100
  real(dp), parameter :: scba_tol = 1.0e-9_dp

  real(dp) :: e_grid(n_en), w_quad(n_en)
  integer :: i

  do i = 1, n_en
    e_grid(i) = e0 + real(i - 1, dp) * de
    w_quad(i) = de
  end do

  ! =========================================================================
  ! Test 6a: digital current step at defective heterointerface (pl=3, 6 PLs)
  ! Also tests the spectral-mode of get_band_resolved_pl_currents.
  ! =========================================================================
  block
    integer, parameter :: n_pl  = 6
    integer, parameter :: n_bnd = n_pl - 1  ! 5 boundaries

    type(TDefectEnsemble) :: ensemble

    complex(dp) :: g_less(n_orb, n_en)
    complex(dp) :: g_great(n_orb, n_en)
    complex(dp) :: g_ret(n_orb, n_en)
    complex(dp) :: g_ret_0(n_orb, n_en)
    complex(dp) :: sigma_0_less(n_orb, n_en)
    complex(dp) :: sigma_0_great(n_orb, n_en)
    complex(dp) :: sigma_less(n_orb, n_en)
    complex(dp) :: sigma_great(n_orb, n_en)
    complex(dp) :: sigma_ret(n_orb, n_en)
    complex(dp) :: sigma_less_prev(n_orb, n_en)

    real(dp) :: j_el_pl(n_bnd), j_hl_pl(n_bnd), j_tot_pl(n_bnd)
    real(dp) :: j_el_sp(n_bnd), j_hl_sp(n_bnd), j_tot_sp(n_bnd)
    real(dp) :: j_spec_pl(n_bnd, n_en)
    real(dp) :: e_split(n_bnd)
    real(dp) :: rate_net_el_ref
    real(dp) :: f_vb, f_cb, abs_g_ret_sq, f_occ_prev, scba_err
    real(dp) :: tg_less_m, tg_great_m, coll_i
    integer :: iter, m_orb, l_bnd

    logical :: converged

    ! --- Build ensemble: 1 defect at pl=3 ---
    call init_defect_ensemble(ensemble, 1, n_pl)
    associate(spec => ensemble%defects(1))
      call init_defect_spec(spec, 2, [1, 2], [0.5_dp, 0.5_dp], 1)
      spec%eps_d          = eps_d
      spec%kbt            = kbt
      spec%density_weight = 1.0_dp
      spec%pl_index       = 3
      spec%mode           = DEFECT_MODE_DELTA
      spec%ec_cut         = ec_cut_val
      spec%ev_cut         = ev_cut_val
      spec%e_shift(1)    = e_shift
      spec%m_dc_capt(1)  = m_dc_capt_val
      spec%m_dc_em(1)    = m_dc_capt_val * exp(-e_shift / kbt)
      spec%m_dv_em(1)    = m_dv_em_val
      spec%m_dv_capt(1)  = m_dv_em_val * exp(-e_shift / kbt)
    end associate

    ! --- Bare contact Green's functions ---
    do i = 1, n_en
      f_vb = 1.0_dp / (1.0_dp + exp((e_grid(i) - mu_v) / kbt))
      f_cb = 1.0_dp / (1.0_dp + exp((e_grid(i) - mu_c) / kbt))
      g_ret_0(1, i) = 1.0_dp / cmplx(0.0_dp, 0.5_dp * gamma_0, dp)
      g_ret_0(2, i) = 1.0_dp / cmplx(0.0_dp, 0.5_dp * gamma_0, dp)
      sigma_0_less(1, i)  = cmplx(0.0_dp,  gamma_0 * f_vb,            dp)
      sigma_0_great(1, i) = cmplx(0.0_dp, -gamma_0 * (1.0_dp - f_vb), dp)
      sigma_0_less(2, i)  = cmplx(0.0_dp,  gamma_0 * f_cb,            dp)
      sigma_0_great(2, i) = cmplx(0.0_dp, -gamma_0 * (1.0_dp - f_cb), dp)
    end do

    sigma_less(:, :)  = (0.0_dp, 0.0_dp)
    sigma_great(:, :) = (0.0_dp, 0.0_dp)
    sigma_ret(:, :)   = (0.0_dp, 0.0_dp)

    ! --- Coupled SCBA loop ---
    f_occ_prev = 0.5_dp
    converged  = .false.
    scba_err   = 1.0_dp
    iter       = 0

    do while (.not. converged .and. iter < max_scba)
      do i = 1, n_en
        do m_orb = 1, n_orb
          g_ret(m_orb, i) = 1.0_dp / &
              & (1.0_dp / g_ret_0(m_orb, i) - sigma_ret(m_orb, i))
          abs_g_ret_sq = real(g_ret(m_orb, i))**2 &
              & + aimag(g_ret(m_orb, i))**2
          g_less(m_orb, i)  = abs_g_ret_sq * &
              & (sigma_0_less(m_orb, i)  + sigma_less(m_orb, i))
          g_great(m_orb, i) = abs_g_ret_sq * &
              & (sigma_0_great(m_orb, i) + sigma_great(m_orb, i))
        end do
      end do

      call solve_defect_occupation(ensemble%defects(1), e_grid, w_quad, &
          & g_less, g_great)

      sigma_less_prev(:, :) = sigma_less(:, :)

      sigma_less(:, :)  = (0.0_dp, 0.0_dp)
      sigma_great(:, :) = (0.0_dp, 0.0_dp)
      sigma_ret(:, :)   = (0.0_dp, 0.0_dp)

      ensemble%pl_recomb_rate(3) = &
          & ensemble%defects(1)%density_weight * ensemble%defects(1)%rate_net_el
      ensemble%total_recomb_rate = ensemble%pl_recomb_rate(3)

      call compute_extended_defect_selfenergy(ensemble, e_grid, w_quad, &
          & sigma_less, sigma_great, sigma_ret)

      scba_err = abs(ensemble%defects(1)%f_occ - f_occ_prev)
      if (scba_err < scba_tol .and. &
          & maxval(abs(aimag(sigma_less) - aimag(sigma_less_prev))) &
          & < scba_tol) then
        converged = .true.
      end if
      f_occ_prev = ensemble%defects(1)%f_occ
      iter = iter + 1
    end do

    if (.not. converged) then
      write(*, '(a,i0,a,es10.3)') &
          & 'FAIL 6a-scba: did not converge in ', max_scba, &
          & ' iter; scba_err=', scba_err
      error stop 1
    end if

    ! --- Final dressed GF ---
    do i = 1, n_en
      do m_orb = 1, n_orb
        g_ret(m_orb, i) = 1.0_dp / &
            & (1.0_dp / g_ret_0(m_orb, i) - sigma_ret(m_orb, i))
        abs_g_ret_sq = real(g_ret(m_orb, i))**2 + aimag(g_ret(m_orb, i))**2
        g_less(m_orb, i)  = abs_g_ret_sq * &
            & (sigma_0_less(m_orb, i)  + sigma_less(m_orb, i))
        g_great(m_orb, i) = abs_g_ret_sq * &
            & (sigma_0_great(m_orb, i) + sigma_great(m_orb, i))
      end do
    end do

    ! --- Continuity mode: get_band_resolved_pl_currents without j_spec_pl ---
    e_split(:) = 0.0_dp
    call get_band_resolved_pl_currents(ensemble, e_grid, w_quad, e_split, &
        & j_el_pl, j_hl_pl, j_tot_pl)

    rate_net_el_ref = ensemble%defects(1)%density_weight &
        & * ensemble%defects(1)%rate_net_el

    ! Verify 1a: digital step profile (pl=3 boundary at L=3)
    ! L=1,2: ec_cut(0.15) >= e_split(0) => j_el = R_el; ev_cut(-0.15) > 0 => j_hl=0
    ! L=3,4,5: pl=3 <= L and ev_cut(-0.15) <= 0 => j_hl = R_hl; j_el = 0
    do l_bnd = 1, 2
      if (abs(j_el_pl(l_bnd) - rate_net_el_ref) > 1.0e-14_dp) then
        write(*, '(a,i0,a,es14.6,a,es14.6)') &
            & 'FAIL 6a-step j_el L=', l_bnd, &
            & ' got=', j_el_pl(l_bnd), ' expect=', rate_net_el_ref
        error stop 1
      end if
      if (abs(j_hl_pl(l_bnd)) > 1.0e-14_dp) then
        write(*, '(a,i0,a,es14.6)') &
            & 'FAIL 6a-step j_hl L=', l_bnd, ' got=', j_hl_pl(l_bnd)
        error stop 1
      end if
    end do
    do l_bnd = 3, n_bnd
      if (abs(j_el_pl(l_bnd)) > 1.0e-14_dp) then
        write(*, '(a,i0,a,es14.6)') &
            & 'FAIL 6a-step j_el L=', l_bnd, ' got=', j_el_pl(l_bnd)
        error stop 1
      end if
      if (abs(j_hl_pl(l_bnd) - rate_net_el_ref) > 1.0e-14_dp) then
        write(*, '(a,i0,a,es14.6,a,es14.6)') &
            & 'FAIL 6a-step j_hl L=', l_bnd, &
            & ' got=', j_hl_pl(l_bnd), ' expect=', rate_net_el_ref
        error stop 1
      end if
    end do

    ! Verify 1b: total current spatially constant
    if (maxval(abs(j_tot_pl - j_tot_pl(1))) / &
        & (abs(j_tot_pl(1)) + 1.0e-30_dp) >= 1.0e-10_dp) then
      write(*, '(a,es14.6)') &
          & 'FAIL 6a-cons: J_tot not constant; max dev=', &
          & maxval(abs(j_tot_pl - j_tot_pl(1)))
      error stop 1
    end if

    ! --- Spectral mode: build j_spec_pl from collision integrand ---
    ! For the single defect at pl=3, the collision integrand per energy point
    ! is: coll(i) = sum_m [(-Im Sigma>) * Im(G<) - Im(Sigma<) * (-Im G>)] * w_quad
    ! We store it in j_spec_pl(l_bnd, i) for all boundaries to represent the
    ! current that "flows" through pl=3 — identical for all L since there is
    ! only 1 defect.
    j_spec_pl(:, :) = 0.0_dp
    do i = 1, n_en
      coll_i = 0.0_dp
      do m_orb = 1, n_orb
        tg_less_m  =  aimag(g_less(m_orb, i))
        tg_great_m = -aimag(g_great(m_orb, i))
        coll_i = coll_i + &
            & (-aimag(sigma_great(m_orb, i))) * tg_less_m &
            & - aimag(sigma_less(m_orb, i)) * tg_great_m
      end do
      do l_bnd = 1, 2
        if (e_grid(i) > e_split(l_bnd)) j_spec_pl(l_bnd, i) = coll_i
      end do
      do l_bnd = 3, n_bnd
        if (e_grid(i) <= e_split(l_bnd)) j_spec_pl(l_bnd, i) = -coll_i
      end do
    end do

    call get_band_resolved_pl_currents(ensemble, e_grid, w_quad, e_split, &
        & j_el_sp, j_hl_sp, j_tot_sp, j_spec_pl)

    ! Verify 2: spectral-mode matches continuity-mode on every boundary
    if (maxval(abs(j_el_sp - j_el_pl)) > 1.0e-14_dp .or. &
        & maxval(abs(j_hl_sp - j_hl_pl)) > 1.0e-14_dp .or. &
        & maxval(abs(j_tot_sp - j_tot_pl)) > 1.0e-14_dp) then
      write(*, '(a,es14.6,a,es14.6)') &
          & 'FAIL 6a-spec: spectral j_tot=', j_tot_sp(1), &
          & ' continuity j_tot=', j_tot_pl(1)
      error stop 1
    end if

    call destruct_defect_ensemble(ensemble)
  end block

  ! =========================================================================
  ! Test 6b: field-enhanced subgap tunneling recombination (Franz-Keldysh).
  !
  ! Physical model: band edges at e_vb=-0.25 eV and e_cb=+0.25 eV.
  ! Resonances e_c=+0.20 and e_v=-0.20 lie 50 meV inside the bandgap.
  ! The bare retarded GF is a Lorentzian centred at the band edge with
  ! subgap tail broadening gamma_tail:
  !   g_ret_0(1,i) = 1 / ((E - e_vb_edge) + i*0.5*gamma_tail)
  !   g_ret_0(2,i) = 1 / ((E - e_cb_edge) + i*0.5*gamma_tail)
  ! At the resonance energies the subgap spectral weight A = gamma_tail /
  ! ((E-E_edge)^2 + (0.5*gamma_tail)^2):
  !   Run 1 (gamma_tail=0.005): A ~ 0.995 / (0.0025 + 0.0000063) ~ 388
  !   Run 2 (gamma_tail=0.050): A ~ 50.0  / (0.0025 + 0.000625)  ~ 16000
  ! rate_B >> rate_A confirming field-enhanced subgap recombination.
  ! =========================================================================
  block
    real(dp), parameter :: e_cb_edge  =  0.25_dp
    real(dp), parameter :: e_vb_edge  = -0.25_dp

    real(dp) :: rate_a, rate_b
    integer  :: run

    do run = 1, 2
      block
        integer, parameter :: n_pl = 6
        type(TDefectEnsemble) :: ensemble

        complex(dp) :: g_less(n_orb, n_en)
        complex(dp) :: g_great(n_orb, n_en)
        complex(dp) :: g_ret(n_orb, n_en)
        complex(dp) :: g_ret_0(n_orb, n_en)
        complex(dp) :: sigma_0_less(n_orb, n_en)
        complex(dp) :: sigma_0_great(n_orb, n_en)
        complex(dp) :: sigma_less(n_orb, n_en)
        complex(dp) :: sigma_great(n_orb, n_en)
        complex(dp) :: sigma_ret(n_orb, n_en)
        complex(dp) :: sigma_less_prev(n_orb, n_en)

        real(dp) :: gamma_tail, f_vb, f_cb, abs_g_ret_sq
        real(dp) :: f_occ_prev, scba_err
        integer  :: iter, m_orb
        logical  :: converged

        ! Run 1: zero field — small intrinsic tail broadening
        ! Run 2: finite field — Franz-Keldysh-enhanced subgap tail
        if (run == 1) then
          gamma_tail = 0.005_dp
        else
          gamma_tail = 0.050_dp
        end if

        call init_defect_ensemble(ensemble, 1, n_pl)
        associate(spec => ensemble%defects(1))
          call init_defect_spec(spec, 2, [1, 2], [0.5_dp, 0.5_dp], 1)
          spec%eps_d          = eps_d  ! 0.0 eV, fixed at centre PL
          spec%kbt            = kbt
          spec%density_weight = 1.0_dp
          spec%pl_index       = 3
          spec%mode           = DEFECT_MODE_DELTA
          spec%ec_cut         = ec_cut_val
          spec%ev_cut         = ev_cut_val
          spec%e_shift(1)    = e_shift
          spec%m_dc_capt(1)  = m_dc_capt_val
          spec%m_dc_em(1)    = m_dc_capt_val * exp(-e_shift / kbt)
          spec%m_dv_em(1)    = m_dv_em_val
          spec%m_dv_capt(1)  = m_dv_em_val * exp(-e_shift / kbt)
        end associate

        ! Lorentzian DOS centred at band edges with field-dependent tail width
        do i = 1, n_en
          f_vb = 1.0_dp / (1.0_dp + exp((e_grid(i) - mu_v) / kbt))
          f_cb = 1.0_dp / (1.0_dp + exp((e_grid(i) - mu_c) / kbt))
          ! VB: Lorentzian at e_vb_edge = -0.25
          g_ret_0(1, i) = 1.0_dp / cmplx( &
              & e_grid(i) - e_vb_edge, 0.5_dp * gamma_tail, dp)
          ! CB: Lorentzian at e_cb_edge = +0.25
          g_ret_0(2, i) = 1.0_dp / cmplx( &
              & e_grid(i) - e_cb_edge, 0.5_dp * gamma_tail, dp)
          sigma_0_less(1, i)  = cmplx(0.0_dp,  gamma_tail * f_vb,            dp)
          sigma_0_great(1, i) = cmplx(0.0_dp, -gamma_tail * (1.0_dp - f_vb), dp)
          sigma_0_less(2, i)  = cmplx(0.0_dp,  gamma_tail * f_cb,            dp)
          sigma_0_great(2, i) = cmplx(0.0_dp, -gamma_tail * (1.0_dp - f_cb), dp)
        end do

        sigma_less(:, :)  = (0.0_dp, 0.0_dp)
        sigma_great(:, :) = (0.0_dp, 0.0_dp)
        sigma_ret(:, :)   = (0.0_dp, 0.0_dp)

        f_occ_prev = 0.5_dp
        converged  = .false.
        scba_err   = 1.0_dp
        iter       = 0

        do while (.not. converged .and. iter < max_scba)
          do i = 1, n_en
            do m_orb = 1, n_orb
              g_ret(m_orb, i) = 1.0_dp / &
                  & (1.0_dp / g_ret_0(m_orb, i) - sigma_ret(m_orb, i))
              abs_g_ret_sq = real(g_ret(m_orb, i))**2 &
                  & + aimag(g_ret(m_orb, i))**2
              g_less(m_orb, i)  = abs_g_ret_sq * &
                  & (sigma_0_less(m_orb, i)  + sigma_less(m_orb, i))
              g_great(m_orb, i) = abs_g_ret_sq * &
                  & (sigma_0_great(m_orb, i) + sigma_great(m_orb, i))
            end do
          end do

          call solve_defect_occupation(ensemble%defects(1), e_grid, w_quad, &
              & g_less, g_great)

          sigma_less_prev(:, :) = sigma_less(:, :)

          sigma_less(:, :)  = (0.0_dp, 0.0_dp)
          sigma_great(:, :) = (0.0_dp, 0.0_dp)
          sigma_ret(:, :)   = (0.0_dp, 0.0_dp)

          ensemble%pl_recomb_rate(3) = &
              & ensemble%defects(1)%density_weight &
              & * ensemble%defects(1)%rate_net_el
          ensemble%total_recomb_rate = ensemble%pl_recomb_rate(3)

          call compute_extended_defect_selfenergy(ensemble, e_grid, w_quad, &
              & sigma_less, sigma_great, sigma_ret)

          scba_err = abs(ensemble%defects(1)%f_occ - f_occ_prev)
          if (scba_err < scba_tol .and. &
              & maxval(abs(aimag(sigma_less) - aimag(sigma_less_prev))) &
              & < scba_tol) then
            converged = .true.
          end if
          f_occ_prev = ensemble%defects(1)%f_occ
          iter = iter + 1
        end do

        if (.not. converged) then
          write(*, '(a,i0,a,i0,a,es10.3)') &
              & 'FAIL 6b-scba run=', run, ' did not converge in ', &
              & max_scba, ' iter; scba_err=', scba_err
          error stop 1
        end if

        if (run == 1) then
          rate_a = ensemble%defects(1)%density_weight &
              & * ensemble%defects(1)%rate_net_el
        else
          rate_b = ensemble%defects(1)%density_weight &
              & * ensemble%defects(1)%rate_net_el
        end if

        call destruct_defect_ensemble(ensemble)
      end block
    end do

    ! Verify: finite field (run 2) enhances subgap recombination > 5x
    if (rate_b <= 5.0_dp * rate_a) then
      write(*, '(a,es14.6,a,es14.6)') &
          & 'FAIL 6b-field: rate_B=', rate_b, &
          & ' should be > 5 * rate_A=', rate_a
      error stop 1
    end if
  end block

  write(*, '(a)') 'Task 6: test_defect_hetero PASSED'

end program test_defect_hetero
