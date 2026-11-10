!!--------------------------------------------------------------------------!
!! libNEGF: a general library for Non-Equilibrium Greens functions.         !
!! Copyright (C) 2012 - 2026                                                !
!!                                                                          !
!! This file is part of libNEGF: a library for                              !
!! Non Equilibrium Green's Functions calculations                           !
!!                                                                          !
!! Developers: Alessandro Pecchia, Daniele Soccodato                        !
!! Former Contributors: Gabriele Penazzi, Luca Latessa, Aldo Di Carlo       !
!!                                                                          !
!! libNEGF is free software: you can redistribute and/or modify it          !
!! under the terms of the GNU Lesser General Public License as published    !
!! by the Free Software Foundation, either version 3 of the License, or     !
!! (at your option) any later version.                                      !
!!                                                                          !
!!  You should have received a copy of the GNU Lesser General Public        !
!!  License along with libNEGF.  If not, see                                !
!!  <http://www.gnu.org/licenses/>.                                         !
!!--------------------------------------------------------------------------!

#:include 'assert.fypp'

!> SRH defect-assisted recombination self-energies for NEGF.
!!
!! Sign convention used throughout:
!!   libnegf stores G< as the pure Keldysh function:
!!     Im(G<_mm) > 0,  Im(G>_mm) < 0.
!!   Tilde quantities (real, non-negative on the diagonal):
!!     tilde_G< = Im(G<) >= 0,  tilde_G> = -Im(G>) >= 0.
!!   Self-energies:
!!     Sigma< = +i * tilde_Sigma<  =>  Im(Sigma<) = +tilde > 0
!!     Sigma> = -i * tilde_Sigma>  =>  Im(Sigma>) = -tilde < 0
!!   Retarded (causality: Im(Sigma^R) <= 0):
!!     Im(Sigma^R) = +0.5*(Im(Sigma>) - Im(Sigma<)) <= 0
module ln_defect_srh

  use ln_precision, only : dp
  use ln_defect_types, only : TDefectSpec, TDefectEnsemble, &
      & DEFECT_MODE_DYSON, init_defect_spec_dyson

  implicit none
  private

  public :: compute_energy_delta_weights
  public :: solve_defect_occupation
  public :: compute_extended_defect_selfenergy
  public :: solve_defect_dyson_keldysh

contains

  !> Compute normalised projection weights delta_w(i) for an energy target.
  !!
  !! Conservation identity: sum_i e_weights(i) * delta_w(i) = 1.0
  !! Out-of-range targets always return all zeros.
  !!
  !! @param[in]  e_grid         Sorted ascending energy grid, size n_en
  !! @param[in]  e_weights      Quadrature weights, size n_en
  !! @param[in]  e_target       Target energy (eV)
  !! @param[in]  eta_smooth     Lorentzian broadening (eV); used only when
  !!                            use_lorentzian=.true.
  !! @param[in]  use_lorentzian Use Lorentzian instead of hat interpolation
  !! @param[out] delta_w        Projection weights, size n_en
  subroutine compute_energy_delta_weights(e_grid, e_weights, e_target, &
      & eta_smooth, use_lorentzian, delta_w)
    real(dp), intent(in) :: e_grid(:)
    real(dp), intent(in) :: e_weights(:)
    real(dp), intent(in) :: e_target
    real(dp), intent(in) :: eta_smooth
    logical, intent(in) :: use_lorentzian
    real(dp), intent(out) :: delta_w(:)

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer :: n_en, j, idx
    real(dp) :: alpha, s, lor_j

    n_en = size(e_grid)
    @:ASSERT(size(e_weights) == n_en)
    @:ASSERT(size(delta_w) == n_en)

    delta_w(:) = 0.0_dp

    ! Out-of-range check before any Lorentzian so tails are never normalised
    if (e_target < e_grid(1) .or. e_target > e_grid(n_en)) return

    if (use_lorentzian) then
      @:ASSERT(eta_smooth > 0.0_dp)
      s = 0.0_dp
      do j = 1, n_en
        lor_j = eta_smooth / &
            & ((e_grid(j) - e_target)**2 + eta_smooth**2) / pi
        delta_w(j) = lor_j
        s = s + e_weights(j) * lor_j
      end do
      if (s > 1.0e-12_dp) then
        do j = 1, n_en
          delta_w(j) = delta_w(j) / s
        end do
        return
      end if
      ! Sum too small: fall through to hat
      delta_w(:) = 0.0_dp
    end if

    ! Hat (linear) interpolation — exact discrete conservation
    idx = 1
    do j = 1, n_en - 1
      if (e_grid(j) <= e_target .and. e_target <= e_grid(j + 1)) then
        idx = j
        exit
      end if
    end do

    alpha = (e_target - e_grid(idx)) / (e_grid(idx + 1) - e_grid(idx))
    delta_w(idx) = (1.0_dp - alpha) / e_weights(idx)
    delta_w(idx + 1) = alpha / e_weights(idx + 1)

  end subroutine compute_energy_delta_weights

  !> Solve for the non-equilibrium defect occupation f_occ using G< and G>
  !! on the energy grid.  Always uses hat interpolation (use_lor = .false.)
  !! so discrete energy conservation holds exactly, including when called
  !! from solve_defect_dyson_keldysh.
  !!
  !! @param[inout] spec         Updated: f_occ, gamma_*, sigma_*, rate_*
  !! @param[in]   e_grid        Energy grid, size n_en
  !! @param[in]   w_quad        Quadrature weights, size n_en
  !! @param[in]   g_less_diag   G<(orb, e), shape (n_orb, n_en)
  !! @param[in]   g_great_diag  G>(orb, e), shape (n_orb, n_en)
  subroutine solve_defect_occupation(spec, e_grid, w_quad, &
      & g_less_diag, g_great_diag)
    type(TDefectSpec), intent(inout) :: spec
    real(dp), intent(in) :: e_grid(:)
    real(dp), intent(in) :: w_quad(:)
    complex(dp), intent(in) :: g_less_diag(:, :)
    complex(dp), intent(in) :: g_great_diag(:, :)

    integer :: n_en, n_orbs, l, i, m_idx, m
    real(dp) :: e_c, e_v
    real(dp) :: tilde_sigma_dc_less, tilde_sigma_dc_great
    real(dp) :: tilde_sigma_dv_less, tilde_sigma_dv_great
    real(dp) :: gamma_dc, gamma_dv, denom, w
    real(dp), allocatable :: delta_w(:)
    real(dp), allocatable :: tg_less(:), tg_great(:)

    n_en = size(e_grid)
    n_orbs = size(spec%orb_indices)
    @:ASSERT(size(w_quad) == n_en)
    @:ASSERT(size(g_less_diag, 2) == n_en)
    @:ASSERT(size(g_great_diag, 2) == n_en)
    @:ASSERT(spec%n_modes > 0)

    tilde_sigma_dc_less  = 0.0_dp
    tilde_sigma_dc_great = 0.0_dp
    tilde_sigma_dv_less  = 0.0_dp
    tilde_sigma_dv_great = 0.0_dp

    allocate(delta_w(n_en))
    allocate(tg_less(n_en))
    allocate(tg_great(n_en))
    tg_less(:)  = 0.0_dp
    tg_great(:) = 0.0_dp
    do m_idx = 1, n_orbs
      m = spec%orb_indices(m_idx)
      w = spec%orb_weights(m_idx)
      do i = 1, n_en
        tg_less(i)  = tg_less(i)  + w * aimag(g_less_diag(m, i))
        tg_great(i) = tg_great(i) - w * aimag(g_great_diag(m, i))
      end do
    end do

    do l = 1, spec%n_modes
      ! Conduction-band resonance
      e_c = spec%eps_d + spec%e_shift(l)
      if (e_c >= spec%ec_cut) then
        call compute_energy_delta_weights(e_grid, w_quad, e_c, &
            & spec%eta_d, .false., delta_w)
        do i = 1, n_en
          tilde_sigma_dc_less  = tilde_sigma_dc_less &
              & + spec%m_dc_capt(l) * w_quad(i) * delta_w(i) * tg_less(i)
          tilde_sigma_dc_great = tilde_sigma_dc_great &
              & + spec%m_dc_em(l) * w_quad(i) * delta_w(i) * tg_great(i)
        end do
      end if

      ! Valence-band resonance
      e_v = spec%eps_d - spec%e_shift(l)
      if (e_v <= spec%ev_cut) then
        call compute_energy_delta_weights(e_grid, w_quad, e_v, &
            & spec%eta_d, .false., delta_w)
        do i = 1, n_en
          tilde_sigma_dv_less  = tilde_sigma_dv_less &
              & + spec%m_dv_capt(l) * w_quad(i) * delta_w(i) * tg_less(i)
          tilde_sigma_dv_great = tilde_sigma_dv_great &
              & + spec%m_dv_em(l) * w_quad(i) * delta_w(i) * tg_great(i)
        end do
      end if
    end do

    deallocate(delta_w)
    deallocate(tg_less)
    deallocate(tg_great)

    gamma_dc = tilde_sigma_dc_less + tilde_sigma_dc_great
    gamma_dv = tilde_sigma_dv_less + tilde_sigma_dv_great
    denom = gamma_dc + gamma_dv
    if (denom < 1.0e-30_dp) then
      spec%f_occ = 0.5_dp
    else
      spec%f_occ = (tilde_sigma_dc_less + tilde_sigma_dv_less) / denom
    end if
    spec%f_occ = min(1.0_dp, max(0.0_dp, spec%f_occ))

    spec%gamma_dc = gamma_dc
    spec%gamma_dv = gamma_dv

    spec%sigma_dc_less  = cmplx(0.0_dp,  tilde_sigma_dc_less,  dp)
    spec%sigma_dc_great = cmplx(0.0_dp, -tilde_sigma_dc_great, dp)
    spec%sigma_dv_less  = cmplx(0.0_dp,  tilde_sigma_dv_less,  dp)
    spec%sigma_dv_great = cmplx(0.0_dp, -tilde_sigma_dv_great, dp)

    spec%rate_net_el = tilde_sigma_dc_less * (1.0_dp - spec%f_occ) &
        & - tilde_sigma_dc_great * spec%f_occ
    spec%rate_net_hl = tilde_sigma_dv_great * spec%f_occ &
        & - tilde_sigma_dv_less * (1.0_dp - spec%f_occ)

  end subroutine solve_defect_occupation

  !> Accumulate SRH defect self-energy onto diagonal arrays.
  !!
  !! Supports both Mode A (DEFECT_MODE_DELTA, uses f_occ and eps_d) and
  !! Mode B (DEFECT_MODE_DYSON, convolves g_d_less/great over omega_grid).
  !!
  !! Output sign convention:
  !!   sigma_less_diag  : Im > 0
  !!   sigma_great_diag : Im < 0
  !!   sigma_ret_diag   : Im <= 0  (causal)
  subroutine compute_extended_defect_selfenergy(ensemble, e_grid, w_quad, &
      & sigma_less_diag, sigma_great_diag, sigma_ret_diag)
    type(TDefectEnsemble), intent(in) :: ensemble
    real(dp), intent(in) :: e_grid(:)
    real(dp), intent(in) :: w_quad(:)
    complex(dp), intent(inout) :: sigma_less_diag(:, :)
    complex(dp), intent(inout) :: sigma_great_diag(:, :)
    complex(dp), intent(inout) :: sigma_ret_diag(:, :)

    integer :: nu, l, i, m_idx, m, n_en, n_orbs_nu, q
    real(dp) :: e_c, e_v, f_d, n_d, wt
    real(dp) :: d_less_im, d_great_im
    real(dp) :: w_q, tg_d_less_q, tg_d_great_q
    real(dp), allocatable :: delta_w(:)

    n_en = size(e_grid)
    allocate(delta_w(n_en))

    do nu = 1, ensemble%num_defects
      associate(spec => ensemble%defects(nu))
        n_d = spec%density_weight
        n_orbs_nu = size(spec%orb_indices)

        if (spec%mode == DEFECT_MODE_DYSON .and. allocated(spec%g_d_less)) then

          ! --- Mode B: convolve g_d_less/great over local omega grid ---
          do q = 1, spec%n_omega
            w_q        = spec%omega_weights(q) / (2.0_dp * acos(-1.0_dp))
            tg_d_less_q  =  aimag(spec%g_d_less(q))   ! tilde_G_d< >= 0
            tg_d_great_q = -aimag(spec%g_d_great(q))  ! tilde_G_d> >= 0

            do l = 1, spec%n_modes

              ! CB resonance from omega_q
              e_c = spec%omega_grid(q) + spec%e_shift(l)
              if (e_c >= spec%ec_cut) then
                call compute_energy_delta_weights(e_grid, w_quad, e_c, &
                    & spec%eta_d, .false., delta_w)
                do i = 1, n_en
                  if (delta_w(i) == 0.0_dp) cycle
                  do m_idx = 1, n_orbs_nu
                    m = spec%orb_indices(m_idx)
                    wt = n_d * spec%orb_weights(m_idx)
                    d_less_im  =  wt * spec%m_dc_em(l) &
                        & * (w_q * tg_d_less_q) * delta_w(i)
                    d_great_im = -wt * spec%m_dc_capt(l) &
                        & * (w_q * tg_d_great_q) * delta_w(i)
                    sigma_less_diag(m, i)  = sigma_less_diag(m, i) &
                        & + cmplx(0.0_dp, d_less_im, dp)
                    sigma_great_diag(m, i) = sigma_great_diag(m, i) &
                        & + cmplx(0.0_dp, d_great_im, dp)
                    sigma_ret_diag(m, i) = sigma_ret_diag(m, i) &
                        & + cmplx(0.0_dp, 0.5_dp*(d_great_im - d_less_im), dp)
                  end do
                end do
              end if

              ! VB resonance from omega_q
              e_v = spec%omega_grid(q) - spec%e_shift(l)
              if (e_v <= spec%ev_cut) then
                call compute_energy_delta_weights(e_grid, w_quad, e_v, &
                    & spec%eta_d, .false., delta_w)
                do i = 1, n_en
                  if (delta_w(i) == 0.0_dp) cycle
                  do m_idx = 1, n_orbs_nu
                    m = spec%orb_indices(m_idx)
                    wt = n_d * spec%orb_weights(m_idx)
                    d_less_im  =  wt * spec%m_dv_em(l) &
                        & * (w_q * tg_d_less_q) * delta_w(i)
                    d_great_im = -wt * spec%m_dv_capt(l) &
                        & * (w_q * tg_d_great_q) * delta_w(i)
                    sigma_less_diag(m, i)  = sigma_less_diag(m, i) &
                        & + cmplx(0.0_dp, d_less_im, dp)
                    sigma_great_diag(m, i) = sigma_great_diag(m, i) &
                        & + cmplx(0.0_dp, d_great_im, dp)
                    sigma_ret_diag(m, i) = sigma_ret_diag(m, i) &
                        & + cmplx(0.0_dp, 0.5_dp*(d_great_im - d_less_im), dp)
                  end do
                end do
              end if

            end do ! modes
          end do ! omega grid

        else

          ! --- Mode A: single-energy delta, use f_occ and eps_d ---
          f_d = spec%f_occ

          do l = 1, spec%n_modes

            e_c = spec%eps_d + spec%e_shift(l)
            if (e_c >= spec%ec_cut) then
              call compute_energy_delta_weights(e_grid, w_quad, e_c, &
                  & spec%eta_d, .false., delta_w)
              do i = 1, n_en
                if (delta_w(i) == 0.0_dp) cycle
                do m_idx = 1, n_orbs_nu
                  m = spec%orb_indices(m_idx)
                  wt = n_d * spec%orb_weights(m_idx)
                  d_less_im  =  wt * spec%m_dc_em(l)   * f_d           * delta_w(i)
                  d_great_im = -wt * spec%m_dc_capt(l) * (1.0_dp - f_d) * delta_w(i)
                  sigma_less_diag(m, i)  = sigma_less_diag(m, i) &
                      & + cmplx(0.0_dp, d_less_im, dp)
                  sigma_great_diag(m, i) = sigma_great_diag(m, i) &
                      & + cmplx(0.0_dp, d_great_im, dp)
                  sigma_ret_diag(m, i) = sigma_ret_diag(m, i) &
                      & + cmplx(0.0_dp, 0.5_dp*(d_great_im - d_less_im), dp)
                end do
              end do
            end if

            e_v = spec%eps_d - spec%e_shift(l)
            if (e_v <= spec%ev_cut) then
              call compute_energy_delta_weights(e_grid, w_quad, e_v, &
                  & spec%eta_d, .false., delta_w)
              do i = 1, n_en
                if (delta_w(i) == 0.0_dp) cycle
                do m_idx = 1, n_orbs_nu
                  m = spec%orb_indices(m_idx)
                  wt = n_d * spec%orb_weights(m_idx)
                  d_less_im  =  wt * spec%m_dv_em(l)   * f_d           * delta_w(i)
                  d_great_im = -wt * spec%m_dv_capt(l) * (1.0_dp - f_d) * delta_w(i)
                  sigma_less_diag(m, i)  = sigma_less_diag(m, i) &
                      & + cmplx(0.0_dp, d_less_im, dp)
                  sigma_great_diag(m, i) = sigma_great_diag(m, i) &
                      & + cmplx(0.0_dp, d_great_im, dp)
                  sigma_ret_diag(m, i) = sigma_ret_diag(m, i) &
                      & + cmplx(0.0_dp, 0.5_dp*(d_great_im - d_less_im), dp)
                end do
              end do
            end if

          end do ! modes

        end if

      end associate
    end do ! defects

    deallocate(delta_w)

  end subroutine compute_extended_defect_selfenergy

  !> Solve the Dyson-Keldysh defect Green's function (Mode B, Aeberhard 2019
  !! Eqs. 10-11).  Builds a local omega grid around eps_d, calls
  !! solve_defect_occupation at each point, constructs the renormalized
  !! G_d^R, G_d^<, G_d^>, and computes spectral-averaged f_occ and net rates.
  !!
  !! Key identity enforced: tilde_G_d<(q) + tilde_G_d>(q) = A_d(q) at every q.
  !!
  !! @param[inout] spec         Updated: g_d_ret/less/great, f_occ, rate_*,
  !!                            sigma_*, gamma_*
  !! @param[in]   e_grid        Extended energy grid, size n_en
  !! @param[in]   w_quad        Quadrature weights, size n_en
  !! @param[in]   g_less_diag   G<(orb, e), shape (n_orb, n_en)
  !! @param[in]   g_great_diag  G>(orb, e), shape (n_orb, n_en)
  subroutine solve_defect_dyson_keldysh(spec, e_grid, w_quad, &
      & g_less_diag, g_great_diag)
    type(TDefectSpec), intent(inout) :: spec
    real(dp), intent(in) :: e_grid(:)
    real(dp), intent(in) :: w_quad(:)
    complex(dp), intent(in) :: g_less_diag(:, :)
    complex(dp), intent(in) :: g_great_diag(:, :)

    real(dp), parameter :: pi = acos(-1.0_dp)

    integer :: q, n_om
    real(dp) :: eps_d_saved, f_occ_saved, d_omega
    real(dp) :: gamma_tot_q, a_d_q, z_norm, w_q
    real(dp) :: tilde_tot_less_q, tilde_tot_great_q
    real(dp) :: tg_d_less_q, tg_d_great_q, f_d_q
    real(dp) :: integral_less, rate_el, rate_hl
    real(dp) :: avg_dc_less, avg_dc_great, avg_dv_less, avg_dv_great

    real(dp), allocatable :: tilde_dc_less(:), tilde_dc_great(:)
    real(dp), allocatable :: tilde_dv_less(:), tilde_dv_great(:)

    n_om = spec%n_omega
    @:ASSERT(n_om >= 2)

    ! Auto-allocate omega arrays if needed
    if (.not. allocated(spec%omega_grid)) then
      call init_defect_spec_dyson(spec)
    else if (size(spec%omega_grid) /= n_om) then
      call init_defect_spec_dyson(spec)
    end if

    allocate(tilde_dc_less(n_om))
    allocate(tilde_dc_great(n_om))
    allocate(tilde_dv_less(n_om))
    allocate(tilde_dv_great(n_om))

    eps_d_saved = spec%eps_d
    f_occ_saved = spec%f_occ

    ! --- Step 1: Build local omega grid with trapezoidal weights ---
    d_omega = 2.0_dp * spec%omega_halfwidth / real(n_om - 1, dp)
    do q = 1, n_om
      spec%omega_grid(q) = eps_d_saved - spec%omega_halfwidth &
          & + real(q - 1, dp) * d_omega
      if (q == 1 .or. q == n_om) then
        spec%omega_weights(q) = 0.5_dp * d_omega
      else
        spec%omega_weights(q) = d_omega
      end if
    end do

    ! --- Pass 1: collect tilde self-energies and build G_d^R ---
    do q = 1, n_om
      spec%eps_d = spec%omega_grid(q)
      call solve_defect_occupation(spec, e_grid, w_quad, &
          & g_less_diag, g_great_diag)
      ! tilde_Sigma< = +Im(Sigma<) >= 0,  tilde_Sigma> = -Im(Sigma>) >= 0
      tilde_dc_less(q)  =  aimag(spec%sigma_dc_less)
      tilde_dc_great(q) = -aimag(spec%sigma_dc_great)
      tilde_dv_less(q)  =  aimag(spec%sigma_dv_less)
      tilde_dv_great(q) = -aimag(spec%sigma_dv_great)

      gamma_tot_q = tilde_dc_less(q) + tilde_dc_great(q) &
          & + tilde_dv_less(q) + tilde_dv_great(q)

      ! G_d^R(q) = 1 / ((omega_q - eps_d) + i*(eta_d + 0.5*gamma_tot))
      ! Im(denom) = eta_d + 0.5*gamma_tot >= eta_d > 0 always (causal)
      spec%g_d_ret(q) = 1.0_dp / cmplx( &
          & spec%omega_grid(q) - eps_d_saved, &
          & spec%eta_d + 0.5_dp * gamma_tot_q, dp)
    end do
    spec%eps_d = eps_d_saved

    ! --- Step 2: Normalize G_d^R (sum rule int A_d dω/2π = 1) ---
    z_norm = 0.0_dp
    do q = 1, n_om
      a_d_q = -2.0_dp * aimag(spec%g_d_ret(q))
      z_norm = z_norm + (spec%omega_weights(q) / (2.0_dp * pi)) * a_d_q
    end do
    @:ASSERT(z_norm > 1.0e-30_dp)

    do q = 1, n_om
      spec%g_d_ret(q) = spec%g_d_ret(q) / z_norm
    end do

    ! --- Pass 2: Keldysh propagators using local occupation f_d(q) ---
    ! Identity enforced: tilde_G_d<(q) + tilde_G_d>(q) = A_d(q) at every q.
    integral_less = 0.0_dp
    rate_el       = 0.0_dp
    rate_hl       = 0.0_dp
    avg_dc_less   = 0.0_dp
    avg_dc_great  = 0.0_dp
    avg_dv_less   = 0.0_dp
    avg_dv_great  = 0.0_dp

    do q = 1, n_om
      w_q   = spec%omega_weights(q) / (2.0_dp * pi)
      a_d_q = -2.0_dp * aimag(spec%g_d_ret(q))   ! normalized, >= 0

      tilde_tot_less_q  = tilde_dc_less(q)  + tilde_dv_less(q)
      tilde_tot_great_q = tilde_dc_great(q) + tilde_dv_great(q)
      gamma_tot_q       = tilde_tot_less_q  + tilde_tot_great_q

      if (gamma_tot_q > 1.0e-30_dp) then
        f_d_q = tilde_tot_less_q / gamma_tot_q
      else
        f_d_q = f_occ_saved
      end if
      f_d_q = min(1.0_dp, max(0.0_dp, f_d_q))

      ! Keldysh: tilde_G_d< = f_d * A_d,  tilde_G_d> = (1-f_d) * A_d
      tg_d_less_q  = f_d_q * a_d_q
      tg_d_great_q = (1.0_dp - f_d_q) * a_d_q

      ! Store in libnegf convention: G_d< = i*tilde, G_d> = -i*tilde
      spec%g_d_less(q)  = cmplx(0.0_dp,  tg_d_less_q,  dp)
      spec%g_d_great(q) = cmplx(0.0_dp, -tg_d_great_q, dp)

      ! Accumulate f_occ = int(tilde_G_d<) dω/2π  (equals int f_d*A_d dω/2π)
      integral_less = integral_less + w_q * tg_d_less_q

      ! Net rates (Eqs. 55 & 58 of Aeberhard, integrated over omega):
      rate_el = rate_el + w_q * ( &
          & tilde_dc_less(q) * tg_d_great_q - tilde_dc_great(q) * tg_d_less_q)
      rate_hl = rate_hl + w_q * ( &
          & tilde_dv_great(q) * tg_d_less_q - tilde_dv_less(q) * tg_d_great_q)

      ! Spectral-averaged self-energies for storage
      avg_dc_less  = avg_dc_less  + w_q * a_d_q * tilde_dc_less(q)
      avg_dc_great = avg_dc_great + w_q * a_d_q * tilde_dc_great(q)
      avg_dv_less  = avg_dv_less  + w_q * a_d_q * tilde_dv_less(q)
      avg_dv_great = avg_dv_great + w_q * a_d_q * tilde_dv_great(q)
    end do

    ! f_occ = int tilde_G_d< dω/2π  (sum rule ensures this is in [0,1])
    spec%f_occ = min(1.0_dp, max(0.0_dp, integral_less))
    spec%rate_net_el = rate_el
    spec%rate_net_hl = rate_hl

    ! Store spectral-averaged self-energies (overwrite per-point garbage)
    spec%sigma_dc_less  = cmplx(0.0_dp,  avg_dc_less,  dp)
    spec%sigma_dc_great = cmplx(0.0_dp, -avg_dc_great, dp)
    spec%sigma_dv_less  = cmplx(0.0_dp,  avg_dv_less,  dp)
    spec%sigma_dv_great = cmplx(0.0_dp, -avg_dv_great, dp)
    spec%gamma_dc = avg_dc_less + avg_dc_great
    spec%gamma_dv = avg_dv_less + avg_dv_great

    deallocate(tilde_dc_less)
    deallocate(tilde_dc_great)
    deallocate(tilde_dv_less)
    deallocate(tilde_dv_great)

  end subroutine solve_defect_dyson_keldysh

end module ln_defect_srh
