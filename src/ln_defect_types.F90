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

!> Data structures for SRH defect-assisted recombination in NEGF.
!! Contains configuration parameters, derived types for individual defect
!! specifications and defect ensembles, plus their allocation lifecycle
!! routines.
module ln_defect_types

  use ln_precision, only : dp

  implicit none
  private

  public :: DEFECT_MODE_DELTA, DEFECT_MODE_DYSON
  public :: DEFECT_KERNEL_CUSTOM, DEFECT_KERNEL_MPR
  public :: TDefectSpec, TDefectEnsemble
  public :: init_defect_ensemble, init_defect_spec, init_defect_spec_dyson
  public :: destruct_defect_spec, destruct_defect_ensemble

  !> Defect mode: delta approximation (no spatial extent)
  integer, parameter :: DEFECT_MODE_DELTA = 0
  !> Defect mode: full Dyson embedding
  integer, parameter :: DEFECT_MODE_DYSON = 1

  !> Scattering kernel: custom (user-supplied)
  integer, parameter :: DEFECT_KERNEL_CUSTOM = 0
  !> Scattering kernel: multiphonon recombination (MPR)
  integer, parameter :: DEFECT_KERNEL_MPR = 1

  !> Specification of a single defect centre for SRH-NEGF calculations.
  !! Holds physical parameters, orbital mapping, and result accumulators.
  type :: TDefectSpec
    !> Defect energy level (eV)
    real(dp) :: eps_d = 0.0_dp
    !> Weight factor for this defect in the ensemble density
    real(dp) :: density_weight = 1.0_dp
    !> Thermal energy k_B T (eV)
    real(dp) :: kbt = 0.02585_dp
    !> Conduction-band cut-off energy (eV)
    real(dp) :: ec_cut = 0.0_dp
    !> Valence-band cut-off energy (eV)
    real(dp) :: ev_cut = 0.0_dp
    !> Broadening parameter for the defect Green's function (eV)
    real(dp) :: eta_d = 1.0e-3_dp
    !> Defect interaction mode (DEFECT_MODE_DELTA or DEFECT_MODE_DYSON)
    integer :: mode = DEFECT_MODE_DELTA
    !> Scattering kernel type (DEFECT_KERNEL_CUSTOM or DEFECT_KERNEL_MPR)
    integer :: kernel_type = DEFECT_KERNEL_MPR
    !> Principal-layer index hosting the defect
    integer :: pl_index = 1
    !> Orbital indices in the device basis coupled to this defect
    integer, allocatable :: orb_indices(:)
    !> Coupling weights for each orbital
    real(dp), allocatable :: orb_weights(:)
    !> Number of vibronic modes in the MPR kernel
    integer :: n_modes = 0
    !> Zero-phonon energy (eV)
    real(dp) :: hbar_omega0 = 0.036_dp
    !> Huang-Rhys factor for electron capture
    real(dp) :: huang_rhys_n = 1.0_dp
    !> Huang-Rhys factor for hole capture
    real(dp) :: huang_rhys_p = 1.0_dp
    !> Ratio of radiative to phonon emission rate for electrons
    real(dp) :: rf_rph_ratio_n = 0.0_dp
    !> Ratio of radiative to phonon emission rate for holes
    real(dp) :: rf_rph_ratio_p = 0.0_dp
    !> Zero-phonon conduction-band coupling matrix element (eV)
    real(dp) :: m0_dc = 0.0_dp
    !> Zero-phonon valence-band coupling matrix element (eV)
    real(dp) :: m0_dv = 0.0_dp
    !> Energy shift for each vibronic mode (eV), size n_modes
    real(dp), allocatable :: e_shift(:)
    !> Conduction-band capture matrix elements per mode (eV), size n_modes
    real(dp), allocatable :: m_dc_capt(:)
    !> Conduction-band emission matrix elements per mode (eV), size n_modes
    real(dp), allocatable :: m_dc_em(:)
    !> Valence-band emission matrix elements per mode (eV), size n_modes
    real(dp), allocatable :: m_dv_em(:)
    !> Valence-band capture matrix elements per mode (eV), size n_modes
    real(dp), allocatable :: m_dv_capt(:)
    !> Defect occupancy (0 = empty, 1 = filled, 0.5 = half-filled)
    real(dp) :: f_occ = 0.5_dp
    !> Lesser self-energy from the conduction band (eV)
    complex(dp) :: sigma_dc_less = (0.0_dp, 0.0_dp)
    !> Greater self-energy from the conduction band (eV)
    complex(dp) :: sigma_dc_great = (0.0_dp, 0.0_dp)
    !> Lesser self-energy from the valence band (eV)
    complex(dp) :: sigma_dv_less = (0.0_dp, 0.0_dp)
    !> Greater self-energy from the valence band (eV)
    complex(dp) :: sigma_dv_great = (0.0_dp, 0.0_dp)
    !> Linewidth from conduction-band coupling (eV)
    real(dp) :: gamma_dc = 0.0_dp
    !> Linewidth from valence-band coupling (eV)
    real(dp) :: gamma_dv = 0.0_dp
    !> Net electron recombination rate (1/s)
    real(dp) :: rate_net_el = 0.0_dp
    !> Net hole recombination rate (1/s)
    real(dp) :: rate_net_hl = 0.0_dp
    ! Mode B (Dyson-Keldysh) local defect mesh — allocated by init_defect_spec_dyson
    !> Number of local omega grid points for Dyson-Keldysh integration
    integer :: n_omega = 41
    !> Half-width of the local omega window around eps_d (eV)
    real(dp) :: omega_halfwidth = 0.05_dp
    !> Local omega grid points (eV), size n_omega
    real(dp), allocatable :: omega_grid(:)
    !> Local omega quadrature weights (eV), size n_omega
    real(dp), allocatable :: omega_weights(:)
    !> Retarded defect Green's function on local grid, size n_omega
    complex(dp), allocatable :: g_d_ret(:)
    !> Lesser defect Green's function on local grid, size n_omega
    complex(dp), allocatable :: g_d_less(:)
    !> Greater defect Green's function on local grid, size n_omega
    complex(dp), allocatable :: g_d_great(:)
  end type TDefectSpec

  !> Collection of defect centres together with ensemble-level accumulators.
  type :: TDefectEnsemble
    !> True when the ensemble has been initialised and is ready for use
    logical :: is_active = .false.
    !> Number of defect centres in the ensemble
    integer :: num_defects = 0
    !> Per-defect specifications, size num_defects
    type(TDefectSpec), allocatable :: defects(:)
    !> Radiative recombination rate per principal layer (1/s), size num_pls
    real(dp), allocatable :: pl_recomb_rate(:)
    !> Sum of pl_recomb_rate over all layers
    real(dp) :: total_recomb_rate = 0.0_dp
  end type TDefectEnsemble

contains

  !> Initialise a TDefectEnsemble with num_defects defect slots and
  !! num_pls principal-layer accumulators.
  !!
  !! @param[inout] ensemble   Ensemble to initialise
  !! @param[in]   num_defects Number of defect centres to allocate
  !! @param[in]   num_pls     Number of principal layers
  subroutine init_defect_ensemble(ensemble, num_defects, num_pls)
    type(TDefectEnsemble), intent(inout) :: ensemble
    integer, intent(in) :: num_defects
    integer, intent(in) :: num_pls

    @:ASSERT(num_defects >= 0)
    @:ASSERT(num_pls >= 0)

    allocate(ensemble%defects(num_defects))
    allocate(ensemble%pl_recomb_rate(num_pls))
    ensemble%num_defects = num_defects
    ensemble%is_active = .true.
    ensemble%pl_recomb_rate(:) = 0.0_dp

  end subroutine init_defect_ensemble

  !> Initialise the orbital mapping and vibronic arrays of a TDefectSpec.
  !!
  !! @param[inout] spec         Defect specification to initialise
  !! @param[in]   n_orbs        Number of coupled orbitals
  !! @param[in]   orb_indices   Orbital indices in the device basis, size n_orbs
  !! @param[in]   orb_weights   Coupling weights per orbital, size n_orbs
  !! @param[in]   n_modes       Number of vibronic modes
  subroutine init_defect_spec(spec, n_orbs, orb_indices, orb_weights, n_modes)
    type(TDefectSpec), intent(inout) :: spec
    integer, intent(in) :: n_orbs
    integer, intent(in) :: orb_indices(:)
    real(dp), intent(in) :: orb_weights(:)
    integer, intent(in) :: n_modes

    @:ASSERT(n_orbs >= 0)
    @:ASSERT(size(orb_indices) == n_orbs)
    @:ASSERT(size(orb_weights) == n_orbs)
    @:ASSERT(n_modes >= 0)

    allocate(spec%orb_indices(n_orbs))
    allocate(spec%orb_weights(n_orbs))
    spec%orb_indices(:) = orb_indices
    spec%orb_weights(:) = orb_weights

    allocate(spec%e_shift(n_modes))
    allocate(spec%m_dc_capt(n_modes))
    allocate(spec%m_dc_em(n_modes))
    allocate(spec%m_dv_em(n_modes))
    allocate(spec%m_dv_capt(n_modes))
    spec%n_modes = n_modes

  end subroutine init_defect_spec

  !> Deallocate all allocatable components of a TDefectSpec.
  !!
  !! @param[inout] spec  Defect specification to clean up
  subroutine destruct_defect_spec(spec)
    type(TDefectSpec), intent(inout) :: spec

    if (allocated(spec%orb_indices)) deallocate(spec%orb_indices)
    if (allocated(spec%orb_weights)) deallocate(spec%orb_weights)
    if (allocated(spec%e_shift)) deallocate(spec%e_shift)
    if (allocated(spec%m_dc_capt)) deallocate(spec%m_dc_capt)
    if (allocated(spec%m_dc_em)) deallocate(spec%m_dc_em)
    if (allocated(spec%m_dv_em)) deallocate(spec%m_dv_em)
    if (allocated(spec%m_dv_capt)) deallocate(spec%m_dv_capt)
    if (allocated(spec%omega_grid))    deallocate(spec%omega_grid)
    if (allocated(spec%omega_weights)) deallocate(spec%omega_weights)
    if (allocated(spec%g_d_ret))       deallocate(spec%g_d_ret)
    if (allocated(spec%g_d_less))      deallocate(spec%g_d_less)
    if (allocated(spec%g_d_great))     deallocate(spec%g_d_great)

  end subroutine destruct_defect_spec

  !> Allocate the Mode B (Dyson-Keldysh) local omega grid arrays.
  !! Safe to call multiple times; performs a re-init deallocating existing
  !! allocations first.
  !!
  !! @param[inout] spec  Defect specification with n_omega and omega_halfwidth set
  subroutine init_defect_spec_dyson(spec)
    type(TDefectSpec), intent(inout) :: spec

    if (allocated(spec%omega_grid))    deallocate(spec%omega_grid)
    if (allocated(spec%omega_weights)) deallocate(spec%omega_weights)
    if (allocated(spec%g_d_ret))       deallocate(spec%g_d_ret)
    if (allocated(spec%g_d_less))      deallocate(spec%g_d_less)
    if (allocated(spec%g_d_great))     deallocate(spec%g_d_great)

    allocate(spec%omega_grid(spec%n_omega))
    allocate(spec%omega_weights(spec%n_omega))
    allocate(spec%g_d_ret(spec%n_omega))
    allocate(spec%g_d_less(spec%n_omega))
    allocate(spec%g_d_great(spec%n_omega))

  end subroutine init_defect_spec_dyson

  !> Destroy a TDefectEnsemble: clean up each defect, then free ensemble arrays.
  !!
  !! @param[inout] ensemble  Ensemble to destroy
  subroutine destruct_defect_ensemble(ensemble)
    type(TDefectEnsemble), intent(inout) :: ensemble

    integer :: ii

    if (allocated(ensemble%defects)) then
      do ii = 1, ensemble%num_defects
        call destruct_defect_spec(ensemble%defects(ii))
      end do
      deallocate(ensemble%defects)
    end if

    if (allocated(ensemble%pl_recomb_rate)) deallocate(ensemble%pl_recomb_rate)

    ensemble%is_active = .false.
    ensemble%num_defects = 0

  end subroutine destruct_defect_ensemble

end module ln_defect_types
