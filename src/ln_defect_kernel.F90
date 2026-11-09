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

!> Kernel routines for SRH defect-assisted recombination in NEGF.
!! Provides the scaled modified Bessel function, the Huang-Rhys thermal
!! weight, and routines to populate TDefectSpec with MPR or custom kernels.
module ln_defect_kernel

  use ln_precision, only : dp
  use ln_defect_types, only : TDefectSpec, DEFECT_KERNEL_CUSTOM

  implicit none
  private

  public :: bessel_i_scaled
  public :: huang_rhys_weight
  public :: build_mpr_kernel
  public :: set_custom_defect_kernel

contains

  !> Compute the exponentially scaled modified Bessel function exp(-z)*I_l(z).
  !! Uses Miller's downward recurrence with dynamic starting order.
  !!
  !! @param[in] l   Order (non-negative integer)
  !! @param[in] z   Argument (real, >= 0)
  !! @result    val exp(-z)*I_l(z)
  pure function bessel_i_scaled(l, z) result(val)
    integer, intent(in) :: l
    real(dp), intent(in) :: z
    real(dp) :: val

    integer :: m_start, k
    real(dp), allocatable :: p(:)
    real(dp) :: s_norm

    if (z == 0.0_dp) then
      if (l == 0) then
        val = 1.0_dp
      else
        val = 0.0_dp
      end if
      return
    end if

    m_start = l + max(40, int(2.0_dp * sqrt(40.0_dp * real(l, dp) + z * z)) + 20)

    allocate(p(0:m_start + 1))
    p(m_start + 1) = 0.0_dp
    p(m_start) = 1.0_dp

    do k = m_start, 1, -1
      p(k - 1) = (2.0_dp * real(k, dp) / z) * p(k) + p(k + 1)
      if (abs(p(k - 1)) > 1.0e100_dp) then
        p(0:m_start + 1) = p(0:m_start + 1) * 1.0e-100_dp
      end if
    end do

    s_norm = p(0) + 2.0_dp * sum(p(1:m_start))
    val = p(l) / s_norm

  end function bessel_i_scaled

  !> Compute the Huang-Rhys thermal weight L_s(l) for vibronic mode l.
  !! Works for negative l (emission) via detailed balance.
  !!
  !! @param[in] l           Vibronic mode index (may be negative)
  !! @param[in] s_hr        Huang-Rhys factor
  !! @param[in] hbar_omega0 Phonon energy (eV)
  !! @param[in] kbt         Thermal energy k_B T (eV)
  !! @result    w_l         Thermal weight
  pure function huang_rhys_weight(l, s_hr, hbar_omega0, kbt) result(w_l)
    integer, intent(in) :: l
    real(dp), intent(in) :: s_hr
    real(dp), intent(in) :: hbar_omega0
    real(dp), intent(in) :: kbt
    real(dp) :: w_l

    real(dp) :: f_b, z, tilde_i, log_w

    f_b = 1.0_dp / (exp(hbar_omega0 / kbt) - 1.0_dp)
    z = 2.0_dp * s_hr * sqrt(f_b * (f_b + 1.0_dp))
    tilde_i = bessel_i_scaled(abs(l), z)
    log_w = -s_hr * (sqrt(f_b + 1.0_dp) - sqrt(f_b))**2 &
        & + real(l, dp) * hbar_omega0 / (2.0_dp * kbt)

    if (tilde_i > 0.0_dp) then
      w_l = exp(log_w) * tilde_i
    else
      w_l = 0.0_dp
    end if

  end function huang_rhys_weight

  !> Fill the MPR scattering kernel arrays of a TDefectSpec.
  !! Computes e_shift, m_dc_capt, m_dc_em, m_dv_em, m_dv_capt for each mode.
  !!
  !! @param[inout] spec  Defect specification (n_modes must already be set)
  subroutine build_mpr_kernel(spec)
    type(TDefectSpec), intent(inout) :: spec

    integer :: l
    real(dp) :: l_n, l_p, m_dc, m_dv, e

    do l = 1, spec%n_modes
      e = real(l, dp) * spec%hbar_omega0
      spec%e_shift(l) = e

      l_n = huang_rhys_weight(l, spec%huang_rhys_n, spec%hbar_omega0, spec%kbt)
      l_p = huang_rhys_weight(l, spec%huang_rhys_p, spec%hbar_omega0, spec%kbt)

      m_dc = spec%m0_dc * ( &
          & spec%rf_rph_ratio_n**2 &
          & + (real(l, dp) - spec%huang_rhys_n)**2 / spec%huang_rhys_n &
          & ) * l_n
      spec%m_dc_capt(l) = m_dc
      spec%m_dc_em(l) = m_dc * exp(-e / spec%kbt)

      m_dv = spec%m0_dv * ( &
          & spec%rf_rph_ratio_p**2 &
          & + (real(l, dp) + spec%huang_rhys_p)**2 / spec%huang_rhys_p &
          & ) * l_p
      spec%m_dv_em(l) = m_dv
      spec%m_dv_capt(l) = m_dv * exp(-e / spec%kbt)
    end do

  end subroutine build_mpr_kernel

  !> Set a custom scattering kernel by supplying e_shift, m_dc_capt, and
  !! m_dv_em directly. The reverse arrays are derived via detailed balance.
  !!
  !! @param[inout] spec       Defect specification to update
  !! @param[in]   e_shift     Energy shift per mode (eV), size n_modes
  !! @param[in]   m_dc_capt   Conduction-band capture matrix elements, size n_modes
  !! @param[in]   m_dv_em     Valence-band emission matrix elements, size n_modes
  subroutine set_custom_defect_kernel(spec, e_shift, m_dc_capt, m_dv_em)
    type(TDefectSpec), intent(inout) :: spec
    real(dp), intent(in) :: e_shift(:)
    real(dp), intent(in) :: m_dc_capt(:)
    real(dp), intent(in) :: m_dv_em(:)

    integer :: l
    real(dp) :: e

    @:ASSERT(size(e_shift) == spec%n_modes)
    @:ASSERT(size(m_dc_capt) == spec%n_modes)
    @:ASSERT(size(m_dv_em) == spec%n_modes)

    spec%e_shift(:) = e_shift
    spec%m_dc_capt(:) = m_dc_capt
    spec%m_dv_em(:) = m_dv_em

    do l = 1, spec%n_modes
      e = spec%e_shift(l)
      spec%m_dc_em(l) = spec%m_dc_capt(l) * exp(-e / spec%kbt)
      spec%m_dv_capt(l) = spec%m_dv_em(l) * exp(-e / spec%kbt)
    end do

    spec%kernel_type = DEFECT_KERNEL_CUSTOM

  end subroutine set_custom_defect_kernel

end module ln_defect_kernel
