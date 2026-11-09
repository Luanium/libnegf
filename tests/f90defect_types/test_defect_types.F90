!> Unit test for Task 1: TDefectSpec / TDefectEnsemble allocation lifecycle.
program test_defect_types

  use ln_defect_types
  use ln_precision, only : dp

  implicit none

  type(TDefectEnsemble) :: ensemble

  ! --- initialise ensemble with 2 defects and 4 principal layers ---
  call init_defect_ensemble(ensemble, 2, 4)

  ! --- initialise each defect spec ---
  call init_defect_spec(ensemble%defects(1), 2, [1, 2], [0.5_dp, 0.5_dp], 10)
  call init_defect_spec(ensemble%defects(2), 2, [3, 4], [0.6_dp, 0.4_dp], 10)

  ! --- structural checks ---
  if (ensemble%num_defects /= 2) then
    write(*, '(a)') 'FAIL: ensemble%num_defects /= 2'
    error stop 1
  end if

  if (.not. (ensemble%is_active .eqv. .true.)) then
    write(*, '(a)') 'FAIL: ensemble%is_active is not .true.'
    error stop 1
  end if

  if (size(ensemble%pl_recomb_rate) /= 4) then
    write(*, '(a)') 'FAIL: size(ensemble%pl_recomb_rate) /= 4'
    error stop 1
  end if

  if (size(ensemble%defects(1)%orb_indices) /= 2) then
    write(*, '(a)') 'FAIL: size(defects(1)%orb_indices) /= 2'
    error stop 1
  end if

  if (size(ensemble%defects(1)%m_dc_capt) /= 10) then
    write(*, '(a)') 'FAIL: size(defects(1)%m_dc_capt) /= 10'
    error stop 1
  end if

  if (ensemble%defects(1)%orb_indices(1) /= 1) then
    write(*, '(a)') 'FAIL: defects(1)%orb_indices(1) /= 1'
    error stop 1
  end if

  if (ensemble%defects(1)%orb_indices(2) /= 2) then
    write(*, '(a)') 'FAIL: defects(1)%orb_indices(2) /= 2'
    error stop 1
  end if

  if (abs(ensemble%defects(1)%orb_weights(1) - 0.5_dp) >= 1.0e-15_dp) then
    write(*, '(a)') 'FAIL: defects(1)%orb_weights(1) not equal to 0.5'
    error stop 1
  end if

  if (ensemble%defects(1)%n_modes /= 10) then
    write(*, '(a)') 'FAIL: defects(1)%n_modes /= 10'
    error stop 1
  end if

  ! --- default value checks ---
  if (abs(ensemble%defects(1)%f_occ - 0.5_dp) >= 1.0e-15_dp) then
    write(*, '(a)') 'FAIL: defects(1)%f_occ default /= 0.5'
    error stop 1
  end if

  if (ensemble%defects(1)%mode /= DEFECT_MODE_DELTA) then
    write(*, '(a)') 'FAIL: defects(1)%mode default /= DEFECT_MODE_DELTA (0)'
    error stop 1
  end if

  ! --- destruction ---
  call destruct_defect_ensemble(ensemble)

  if (allocated(ensemble%defects) .eqv. .true.) then
    write(*, '(a)') 'FAIL: ensemble%defects still allocated after destruction'
    error stop 1
  end if

  if (allocated(ensemble%pl_recomb_rate) .eqv. .true.) then
    write(*, '(a)') 'FAIL: ensemble%pl_recomb_rate still allocated after destruction'
    error stop 1
  end if

  write(*, '(a)') 'Task 1: test_defect_types PASSED'

end program test_defect_types
