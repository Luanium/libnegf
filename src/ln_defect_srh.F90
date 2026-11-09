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
!!   Tilde quantities (real, non-negative on the diagonal) are defined as:
!!     tilde_G< = Im(G<) >= 0,  tilde_G> = -Im(G>) >= 0.
!!   Self-energies:
!!     Sigma< = +i * tilde_Sigma<  =>  Im(Sigma<) = +tilde_Sigma< > 0
!!     Sigma> = -i * tilde_Sigma>  =>  Im(Sigma>) = -tilde_Sigma> < 0
!!   Retarded self-energy (causality: Im(Sigma^R) <= 0):
!!     Sigma^R = +1/2 * (Sigma> - Sigma<)
!!     Im(Sigma^R) = +1/2 * (Im(Sigma>) - Im(Sigma<))
!!                 = +1/2 * (-tilde_Sigma> - tilde_Sigma<) <= 0
!!   In terms of d_less_im = +tilde_Sigma< >= 0 and
!!                d_great_im = -tilde_Sigma> <= 0:
!!     Im(Sigma^R) = +0.5 * (d_great_im - d_less_im) <= 0
module ln_defect_srh

  use ln_precision, only : dp
  use ln_defect_types, only : TDefectSpec, TDefectEnsemble, DEFECT_MODE_DYSON

  implicit none
  private

  public :: compute_energy_delta_weights
  public :: solve_defect_occupation
  public :: compute_extended_defect_selfenergy

contains

  !> Compute normalised projection weights delta_w(i) for an energy target.
  !!
  !! The weights satisfy the conservation identity:
  !!   sum_i e_weights(i) * delta_w(i) = 1.0
  !!
  !! Out-of-range targets (e_target outside e_grid) always return all zeros,
  !! regardless of use_lorentzian, to avoid artificially normalising tails.
  !!
  !! @param[in]  e_grid         Energy grid points, size n_en (sorted ascending)
  !! @param[in]  e_weights      Quadrature weights, size n_en
  !! @param[in]  e_target       Target energy (eV)
  !! @param[in]  eta_smooth     Lorentzian broadening (eV); only used when
  !!                            use_lorentzian=.true.
  !! @param[in]  use_lorentzian Use Lorentzian smearing instead of hat
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

    ! Early exit for out-of-range targets (before any Lorentzian computation
    ! so that tails outside the grid are NOT normalised to 1).
    if (e_target < e_grid(1) .or. e_target > e_grid(n_en)) return

    if (use_lorentzian) then
      @:ASSERT(eta_smooth > 0.0_dp)
      ! Lorentzian smearing
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
      ! Lorentzian sum too small: fall through to hat-interpolation
      delta_w(:) = 0.0_dp
    end if

    ! Hat (linear) interpolation — exact conservation by construction
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

  !> Solve for the non-equilibrium defect occupation f_occ using the
  !! NEGF spectral functions G< and G> supplied on the energy grid.
  !!
  !! @param[inout] spec           Defect specification; f_occ, gamma_dc,
  !!                              gamma_dv, sigma_*, rate_* are updated
  !! @param[in]   e_grid          Energy grid, size n_en
  !! @param[in]   w_quad          Quadrature weights, size n_en
  !! @param[in]   g_less_diag     G<(orb, e), shape (n_orb, n_en)
  !! @param[in]   g_great_diag    G>(orb, e), shape (n_orb, n_en)
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
    logical :: use_lor
    real(dp), allocatable :: delta_w(:)
    real(dp), allocatable :: tg_less(:), tg_great(:)

    n_en = size(e_grid)
    n_orbs = size(spec%orb_indices)
    @:ASSERT(size(w_quad) == n_en)
    @:ASSERT(size(g_less_diag, 2) == n_en)
    @:ASSERT(size(g_great_diag, 2) == n_en)
    @:ASSERT(spec%n_modes > 0)

    use_lor = (spec%mode == DEFECT_MODE_DYSON)

    tilde_sigma_dc_less  = 0.0_dp
    tilde_sigma_dc_great = 0.0_dp
    tilde_sigma_dv_less  = 0.0_dp
    tilde_sigma_dv_great = 0.0_dp

    allocate(delta_w(n_en))
    ! Pre-compute orbital-weighted tilde G< and G> for each energy point.
    ! tilde_G< = +Im(G<) >= 0,  tilde_G> = -Im(G>) >= 0
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
      ! --- Conduction-band resonance ---
      e_c = spec%eps_d + spec%e_shift(l)
      if (e_c >= spec%ec_cut) then
        call compute_energy_delta_weights(e_grid, w_quad, e_c, &
            & spec%eta_d, use_lor, delta_w)
        do i = 1, n_en
          tilde_sigma_dc_less  = tilde_sigma_dc_less &
              & + spec%m_dc_capt(l) * w_quad(i) * delta_w(i) * tg_less(i)
          tilde_sigma_dc_great = tilde_sigma_dc_great &
              & + spec%m_dc_em(l) * w_quad(i) * delta_w(i) * tg_great(i)
        end do
      end if

      ! --- Valence-band resonance ---
      e_v = spec%eps_d - spec%e_shift(l)
      if (e_v <= spec%ev_cut) then
        call compute_energy_delta_weights(e_grid, w_quad, e_v, &
            & spec%eta_d, use_lor, delta_w)
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

    ! Defect occupation
    gamma_dc = tilde_sigma_dc_less + tilde_sigma_dc_great
    gamma_dv = tilde_sigma_dv_less + tilde_sigma_dv_great
    denom = gamma_dc + gamma_dv
    if (denom < 1.0e-30_dp) then
      spec%f_occ = 0.5_dp
    else
      spec%f_occ = (tilde_sigma_dc_less + tilde_sigma_dv_less) / denom
    end if
    spec%f_occ = min(1.0_dp, max(0.0_dp, spec%f_occ))

    ! Store linewidths
    spec%gamma_dc = gamma_dc
    spec%gamma_dv = gamma_dv

    ! Store native complex self-energies:
    !   Sigma< = +i * tilde  =>  Im = +tilde > 0
    !   Sigma> = -i * tilde  =>  Im = -tilde < 0
    spec%sigma_dc_less  = cmplx(0.0_dp,  tilde_sigma_dc_less,  dp)
    spec%sigma_dc_great = cmplx(0.0_dp, -tilde_sigma_dc_great, dp)
    spec%sigma_dv_less  = cmplx(0.0_dp,  tilde_sigma_dv_less,  dp)
    spec%sigma_dv_great = cmplx(0.0_dp, -tilde_sigma_dv_great, dp)

    ! Net recombination rates
    spec%rate_net_el = tilde_sigma_dc_less * (1.0_dp - spec%f_occ) &
        & - tilde_sigma_dc_great * spec%f_occ
    spec%rate_net_hl = tilde_sigma_dv_great * spec%f_occ &
        & - tilde_sigma_dv_less * (1.0_dp - spec%f_occ)

  end subroutine solve_defect_occupation

  !> Accumulate SRH defect self-energy contributions onto diagonal arrays.
  !!
  !! Uses the SAME ec_cut / ev_cut guards as solve_defect_occupation to
  !! guarantee exact discrete particle/current conservation.
  !!
  !! Sign convention for output arrays:
  !!   sigma_less_diag  : Im > 0  (Sigma< = +i * tilde_Sigma<)
  !!   sigma_great_diag : Im < 0  (Sigma> = -i * tilde_Sigma>)
  !!   sigma_ret_diag   : Im <= 0 (causality: Im = +0.5*(Im(Sigma>) - Im(Sigma<)))
  !!
  !! @param[in]    ensemble         Defect ensemble
  !! @param[in]    e_grid           Energy grid, size n_en
  !! @param[in]    w_quad           Quadrature weights, size n_en
  !! @param[inout] sigma_less_diag  Lesser self-energy diagonal, (n_orb, n_en)
  !! @param[inout] sigma_great_diag Greater self-energy diagonal, (n_orb, n_en)
  !! @param[inout] sigma_ret_diag   Retarded self-energy diagonal, (n_orb, n_en)
  subroutine compute_extended_defect_selfenergy(ensemble, e_grid, w_quad, &
      & sigma_less_diag, sigma_great_diag, sigma_ret_diag)
    type(TDefectEnsemble), intent(in) :: ensemble
    real(dp), intent(in) :: e_grid(:)
    real(dp), intent(in) :: w_quad(:)
    complex(dp), intent(inout) :: sigma_less_diag(:, :)
    complex(dp), intent(inout) :: sigma_great_diag(:, :)
    complex(dp), intent(inout) :: sigma_ret_diag(:, :)

    integer :: nu, l, i, m_idx, m, n_en, n_orbs_nu
    real(dp) :: e_c, e_v, f_d, n_d, wt
    logical :: use_lor
    real(dp), allocatable :: delta_w(:)
    real(dp) :: d_less_im, d_great_im

    n_en = size(e_grid)
    allocate(delta_w(n_en))

    do nu = 1, ensemble%num_defects
      associate(spec => ensemble%defects(nu))
        f_d = spec%f_occ
        n_d = spec%density_weight
        n_orbs_nu = size(spec%orb_indices)
        use_lor = (spec%mode == DEFECT_MODE_DYSON)

        do l = 1, spec%n_modes

          ! --- Conduction-band resonance (guard matches solve_defect_occupation) ---
          e_c = spec%eps_d + spec%e_shift(l)
          if (e_c >= spec%ec_cut) then
            call compute_energy_delta_weights(e_grid, w_quad, e_c, &
                & spec%eta_d, use_lor, delta_w)
            do i = 1, n_en
              if (delta_w(i) == 0.0_dp) cycle
              do m_idx = 1, n_orbs_nu
                m = spec%orb_indices(m_idx)
                wt = n_d * spec%orb_weights(m_idx)
                ! d_less_im  = +tilde_Sigma< >= 0  (Im(Sigma<) > 0)
                d_less_im  =  wt * spec%m_dc_em(l)   * f_d           * delta_w(i)
                ! d_great_im = -tilde_Sigma> <= 0  (Im(Sigma>) < 0)
                d_great_im = -wt * spec%m_dc_capt(l) * (1.0_dp - f_d) * delta_w(i)
                sigma_less_diag(m, i)  = sigma_less_diag(m, i) &
                    & + cmplx(0.0_dp, d_less_im, dp)
                sigma_great_diag(m, i) = sigma_great_diag(m, i) &
                    & + cmplx(0.0_dp, d_great_im, dp)
                ! Im(Sigma^R) = +0.5*(Im(Sigma>) - Im(Sigma<)) <= 0  (causal)
                sigma_ret_diag(m, i) = sigma_ret_diag(m, i) &
                    & + cmplx(0.0_dp, 0.5_dp * (d_great_im - d_less_im), dp)
              end do
            end do
          end if

          ! --- Valence-band resonance (guard matches solve_defect_occupation) ---
          e_v = spec%eps_d - spec%e_shift(l)
          if (e_v <= spec%ev_cut) then
            call compute_energy_delta_weights(e_grid, w_quad, e_v, &
                & spec%eta_d, use_lor, delta_w)
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
                    & + cmplx(0.0_dp, 0.5_dp * (d_great_im - d_less_im), dp)
              end do
            end do
          end if

        end do ! modes
      end associate
    end do ! defects

    deallocate(delta_w)

  end subroutine compute_extended_defect_selfenergy

end module ln_defect_srh
