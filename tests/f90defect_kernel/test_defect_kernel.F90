!> Unit test for Task 2: ln_defect_kernel routines.
!! Tests bessel_i_scaled, huang_rhys_weight, and build_mpr_kernel.
program test_defect_kernel

  use ln_defect_kernel
  use ln_precision, only : dp
  use ln_defect_types

  implicit none

  ! --- Test 2a: Bessel partition-of-unity and spot value ---
  block
    real(dp) :: z, s
    integer :: l
    real(dp), parameter :: z_vals(3) = [0.5_dp, 5.0_dp, 50.0_dp]
    integer :: iv

    do iv = 1, 3
      z = z_vals(iv)
      s = bessel_i_scaled(0, z)
      do l = 1, 80
        s = s + 2.0_dp * bessel_i_scaled(l, z)
      end do
      if (abs(s - 1.0_dp) >= 1.0e-13_dp) then
        write(*, '(a,es12.4,a,es12.4)') &
            & 'FAIL 2a: Bessel sum != 1 for z=', z, '  err=', abs(s - 1.0_dp)
        error stop 1
      end if
    end do

    ! Spot check: I_1(1) known value
    if (abs(bessel_i_scaled(1, 1.0_dp) - exp(-1.0_dp) * 0.565159103992485_dp) &
        & >= 1.0e-12_dp) then
      write(*, '(a)') 'FAIL 2a: bessel_i_scaled(1,1) spot check failed'
      error stop 1
    end if
  end block

  ! --- Test 2b: Huang-Rhys normalisation and detailed balance ---
  block
    real(dp), parameter :: hbar_omega0 = 0.036_dp
    real(dp), parameter :: kbt = 0.02585_dp
    real(dp), parameter :: s_hr = 2.5_dp
    real(dp) :: s, w_pos, w_neg
    integer :: l

    s = 0.0_dp
    do l = -60, 60
      s = s + huang_rhys_weight(l, s_hr, hbar_omega0, kbt)
    end do
    if (abs(s - 1.0_dp) >= 1.0e-12_dp) then
      write(*, '(a,es12.4)') 'FAIL 2b: Huang-Rhys sum != 1, err=', abs(s - 1.0_dp)
      error stop 1
    end if

    do l = 1, 10
      w_pos = huang_rhys_weight(l, s_hr, hbar_omega0, kbt)
      w_neg = huang_rhys_weight(-l, s_hr, hbar_omega0, kbt)
      if (abs(w_neg - w_pos * exp(-real(l, dp) * hbar_omega0 / kbt)) >= 1.0e-13_dp) then
        write(*, '(a,i0)') 'FAIL 2b: detailed balance violation at l=', l
        error stop 1
      end if
    end do
  end block

  ! --- Test 2c: MPR detailed balance ---
  block
    type(TDefectSpec) :: spec
    integer :: l
    real(dp) :: ratio_dc, ratio_dv, expected

    call init_defect_spec(spec, 1, [1], [1.0_dp], 30)
    spec%hbar_omega0 = 0.036_dp
    spec%kbt = 0.02585_dp
    spec%huang_rhys_n = 2.5_dp
    spec%huang_rhys_p = 2.5_dp
    spec%rf_rph_ratio_n = 0.0_dp
    spec%rf_rph_ratio_p = 0.0_dp
    spec%m0_dc = 1.0_dp
    spec%m0_dv = 1.0_dp
    spec%n_modes = 30

    call build_mpr_kernel(spec)

    do l = 1, 30
      expected = exp(-real(l, dp) * spec%hbar_omega0 / spec%kbt)

      ratio_dc = spec%m_dc_em(l) / spec%m_dc_capt(l)
      if (abs(ratio_dc - expected) >= 1.0e-14_dp) then
        write(*, '(a,i0,a,es12.4)') &
            & 'FAIL 2c: m_dc_em/m_dc_capt wrong at l=', l, '  err=', &
            & abs(ratio_dc - expected)
        error stop 1
      end if

      ratio_dv = spec%m_dv_capt(l) / spec%m_dv_em(l)
      if (abs(ratio_dv - expected) >= 1.0e-14_dp) then
        write(*, '(a,i0,a,es12.4)') &
            & 'FAIL 2c: m_dv_capt/m_dv_em wrong at l=', l, '  err=', &
            & abs(ratio_dv - expected)
        error stop 1
      end if
    end do

    call destruct_defect_spec(spec)
  end block

  write(*, '(a)') 'Task 2: test_defect_kernel PASSED'

end program test_defect_kernel
