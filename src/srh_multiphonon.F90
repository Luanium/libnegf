!!--------------------------------------------------------------------------!
!! libNEGF: a general library for Non-Equilibrium Greens functions.         !
!! Copyright (C) 2012 - 2026                                                !
!!                                                                          !
!! SRH multiphonon defect interaction (Stage 7).                            !
!! Physically inelastic (Aeberhard 2019): TInelastic + energy-indexed       !
!! Sigma cache + multiphonon ladder. Occupation fd diagnostics; fd weights  !
!! emission/absorption in Sigma. Full defect GF Gd still deferred.          !
!!--------------------------------------------------------------------------!

module srh_multiphonon

  use ln_precision, only : dp
  use interactions, only : TInteraction
  use ln_inelastic, only : TInelastic
  use ln_structure, only : TStruct_info
  use mat_def, only : c_dns, z_dns, create, mat_destroy => destroy
  use ln_cache, only : TMatrixCacheMem, TMatLabel

  implicit none
  private

  public :: TSrhMultiphonon
  public :: TSrhMultiphonon_create
  public :: TSrhMultiphonon_init

  real(dp), parameter :: HBAR_EV_S = 6.582119569e-16_dp
  real(dp), parameter :: KB_EV = 8.617333262145e-5_dp
  integer,  parameter :: LMAX_DEFAULT = 4
  real(dp), parameter :: WQ_DEFAULT = 0.036_dp
  real(dp), parameter :: S_DEFAULT  = 3.0_dp

  type, extends(TInelastic) :: TSrhMultiphonon

    real(dp) :: Ed0 = 0.0_dp
    real(dp) :: rho_d = 0.0_dp
    real(dp) :: Cn = 0.0_dp
    real(dp) :: Cp = 0.0_dp
    real(dp) :: temperature = 300.0_dp
    integer :: pl_index = 0
    integer :: orb_index = 0
    real(dp) :: fd = 0.5_dp
    real(dp) :: r_net = 0.0_dp
    real(dp) :: n_loc = 0.0_dp
    real(dp) :: n1 = 0.0_dp
    real(dp) :: V2 = 0.0_dp
    real(dp) :: S_hr = S_DEFAULT
    integer  :: l_max = LMAX_DEFAULT
    real(dp) :: Nq = 0.0_dp
    real(dp) :: dE = 0.0_dp
    integer  :: iEshift = 0
    real(dp), allocatable :: Ml(:)

  contains

    procedure :: add_sigma_r_sp
    procedure :: add_sigma_r_dp
    procedure :: add_sigma_n_sp
    procedure :: add_sigma_n_dp
    procedure :: get_sigma_n_blk
    procedure :: get_sigma_n_mat
    procedure :: set_Gr
    procedure :: set_Gn
    procedure :: compute_Sigma_r
    procedure :: compute_Sigma_n
    procedure :: destroy_Sigma_r
    procedure :: destroy_Sigma_n
    procedure :: destroy
    procedure :: set_EnGrid
    procedure :: update_occupation

  end type TSrhMultiphonon

contains

  subroutine TSrhMultiphonon_create(this)
    class(TInteraction), allocatable :: this
    allocate(TSrhMultiphonon::this)
  end subroutine TSrhMultiphonon_create

  subroutine TSrhMultiphonon_init(this, struct, Ed0, rho_d, Cn, Cp, temperature, &
      & pl_index, orb_index, niter)
    type(TSrhMultiphonon), intent(inout) :: this
    type(TStruct_info), intent(in) :: struct
    real(dp), intent(in) :: Ed0, rho_d, Cn, Cp, temperature
    integer, intent(in) :: pl_index, orb_index, niter

    real(dp) :: gamma_eV, x, fact
    integer :: l

    this%descriptor = "SRH multiphonon defect (Stage 7 fd diagnostic)"
    this%struct = struct
    this%Ed0 = Ed0
    this%rho_d = rho_d
    this%Cn = Cn
    this%Cp = Cp
    this%temperature = temperature
    this%pl_index = pl_index
    this%orb_index = orb_index
    this%scba_niter = niter
    this%fd = 0.5_dp
    this%r_net = 0.0_dp
    this%n_loc = 0.0_dp
    this%n1 = 0.0_dp

    this%wq = WQ_DEFAULT
    this%S_hr = S_DEFAULT
    this%l_max = LMAX_DEFAULT

    x = this%wq / (KB_EV * max(temperature, 1.0_dp))
    if (x > 50.0_dp) then
      this%Nq = 0.0_dp
    else
      this%Nq = 1.0_dp / (exp(x) - 1.0_dp)
    end if

    if (allocated(this%Ml)) deallocate(this%Ml)
    allocate(this%Ml(0:this%l_max))
    fact = 1.0_dp
    this%Ml(0) = exp(-this%S_hr)
    do l = 1, this%l_max
      fact = fact * real(l, dp)
      this%Ml(l) = exp(-this%S_hr) * (this%S_hr ** l) / fact
    end do

    ! Mild test floor: large enough for a visible SCBA signal, small enough
    ! not to dominate the whole device current (0.1 eV was far too strong).
    gamma_eV = HBAR_EV_S * max(Cn, Cp) * max(rho_d, 0.0_dp)
    if (gamma_eV < 5.0e-3_dp) gamma_eV = 5.0e-3_dp
    this%V2 = gamma_eV * gamma_eV

    if (.not. associated(this%sigma_r)) then
      allocate(this%sigma_r, source=TMatrixCacheMem(tagname='SRH_Sigma_r'))
    end if
    if (.not. associated(this%sigma_n)) then
      allocate(this%sigma_n, source=TMatrixCacheMem(tagname='SRH_Sigma_n'))
    end if

  end subroutine TSrhMultiphonon_init

  subroutine set_EnGrid(this, dE)
    class(TSrhMultiphonon), intent(inout) :: this
    real(dp), intent(in) :: dE
    this%dE = abs(dE)
    if (this%dE > 0.0_dp) then
      this%iEshift = max(1, nint(this%wq / this%dE))
    else
      this%iEshift = 0
    end if
  end subroutine set_EnGrid

  logical function site_ok(this)
    class(TSrhMultiphonon), intent(in) :: this
    integer :: nbl, pl_start, pl_end, nrow

    site_ok = .false.
    if (this%pl_index < 1 .or. this%orb_index < 1) return
    nbl = this%struct%num_PLs
    if (this%pl_index > nbl) return
    pl_start = this%struct%mat_PL_start(this%pl_index)
    pl_end   = this%struct%mat_PL_end(this%pl_index)
    nrow = pl_end - pl_start + 1
    if (this%orb_index > nrow) return
    site_ok = .true.
  end function site_ok

  subroutine add_sigma_r_sp(this, esh, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(c_dns), dimension(:,:), intent(inout) :: esh
    integer, intent(in), optional :: en_index, k_index, spin
    type(TMatLabel) :: label
    type(z_dns), pointer :: tmp
    integer :: pl

    if (this%scba_iter .eq. 0) return
    if (.not. site_ok(this)) return
    if (.not. associated(this%sigma_r)) return
    if (.not. present(en_index) .or. .not. present(k_index)) return

    pl = this%pl_index
    label%energy_point = en_index
    label%kpoint = k_index
    label%spin = 1
    label%row_block = pl
    label%col_block = pl

    if (.not. this%sigma_r%is_cached(label)) return
    call this%sigma_r%retrieve_pointer(tmp, label)
    esh(pl, pl)%val = esh(pl, pl)%val - tmp%val

  end subroutine add_sigma_r_sp

  subroutine add_sigma_r_dp(this, esh, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(z_dns), dimension(:,:), intent(inout) :: esh
    integer, intent(in), optional :: en_index, k_index, spin
    type(TMatLabel) :: label
    type(z_dns), pointer :: tmp
    integer :: pl

    if (this%scba_iter .eq. 0) return
    if (.not. site_ok(this)) return
    if (.not. associated(this%sigma_r)) return
    if (.not. present(en_index) .or. .not. present(k_index)) return

    pl = this%pl_index
    label%energy_point = en_index
    label%kpoint = k_index
    label%spin = 1
    label%row_block = pl
    label%col_block = pl

    if (.not. this%sigma_r%is_cached(label)) return
    call this%sigma_r%retrieve_pointer(tmp, label)
    esh(pl, pl)%val = esh(pl, pl)%val - tmp%val

  end subroutine add_sigma_r_dp

  subroutine add_sigma_n_sp(this, sigma, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(c_dns), dimension(:,:), intent(inout) :: sigma
    integer, intent(in), optional :: en_index, k_index, spin
    type(TMatLabel) :: label
    type(z_dns), pointer :: tmp
    integer :: pl

    if (this%scba_iter .eq. 0) return
    if (.not. site_ok(this)) return
    if (.not. associated(this%sigma_n)) return
    if (.not. present(en_index) .or. .not. present(k_index)) return

    pl = this%pl_index
    label%energy_point = en_index
    label%kpoint = k_index
    label%spin = 1
    label%row_block = pl
    label%col_block = pl

    if (.not. this%sigma_n%is_cached(label)) return
    call this%sigma_n%retrieve_pointer(tmp, label)
    sigma(pl, pl)%val = sigma(pl, pl)%val + tmp%val

  end subroutine add_sigma_n_sp

  subroutine add_sigma_n_dp(this, sigma, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(z_dns), dimension(:,:), intent(inout) :: sigma
    integer, intent(in), optional :: en_index, k_index, spin
    type(TMatLabel) :: label
    type(z_dns), pointer :: tmp
    integer :: pl

    if (this%scba_iter .eq. 0) return
    if (.not. site_ok(this)) return
    if (.not. associated(this%sigma_n)) return
    if (.not. present(en_index) .or. .not. present(k_index)) return

    pl = this%pl_index
    label%energy_point = en_index
    label%kpoint = k_index
    label%spin = 1
    label%row_block = pl
    label%col_block = pl

    if (.not. this%sigma_n%is_cached(label)) return
    call this%sigma_n%retrieve_pointer(tmp, label)
    sigma(pl, pl)%val = sigma(pl, pl)%val + tmp%val

  end subroutine add_sigma_n_dp

  subroutine get_sigma_n_blk(this, blk_sigma_n, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(z_dns), dimension(:,:), intent(inout) :: blk_sigma_n
    integer, intent(in), optional :: en_index, k_index, spin
    type(TMatLabel) :: label
    integer :: pl

    if (this%scba_iter .eq. 0) return
    if (.not. site_ok(this)) return
    if (.not. associated(this%sigma_n)) return
    if (.not. present(en_index) .or. .not. present(k_index)) return

    pl = this%pl_index
    label%energy_point = en_index
    label%kpoint = k_index
    label%spin = 1
    label%row_block = pl
    label%col_block = pl
    if (this%sigma_n%is_cached(label)) then
      call this%sigma_n%retrieve(blk_sigma_n(pl, pl), label)
    end if

  end subroutine get_sigma_n_blk

  subroutine get_sigma_n_mat(this, sigma_n, ii, jj, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(z_dns), intent(inout) :: sigma_n
    integer, intent(in) :: ii, jj
    integer, intent(in), optional :: en_index, k_index, spin
    type(TMatLabel) :: label

    if (this%scba_iter .eq. 0) return
    if (.not. site_ok(this)) return
    if (ii /= this%pl_index .or. jj /= this%pl_index) return
    if (.not. associated(this%sigma_n)) return
    if (.not. present(en_index) .or. .not. present(k_index)) return

    label%energy_point = en_index
    label%kpoint = k_index
    label%spin = 1
    label%row_block = ii
    label%col_block = jj
    if (this%sigma_n%is_cached(label)) then
      call this%sigma_n%retrieve(sigma_n, label)
    end if

  end subroutine get_sigma_n_mat

  subroutine set_Gr(this, Gr, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(z_dns), dimension(:,:), intent(in) :: Gr
    integer, intent(in), optional :: en_index, k_index, spin
  end subroutine set_Gr

  subroutine set_Gn(this, Gn, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    type(z_dns), dimension(:,:), intent(in) :: Gn
    integer, intent(in), optional :: en_index, k_index, spin
  end subroutine set_Gn

  subroutine compute_Sigma_r(this, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    integer, intent(in), optional :: en_index, k_index, spin

    integer :: nE, nK, nbl, pl, orb, nrow, iE, iK, pl_start, pl_end, l, iEp, iEm
    type(TMatLabel) :: label, labelp
    type(z_dns), pointer :: pG
    type(z_dns) :: sig
    complex(dp) :: acc
    real(dp) :: fac_em, fac_abs, wgt

    if (.not. site_ok(this)) return
    if (this%V2 <= 0.0_dp) return
    if (.not. associated(this%G_r)) return
    if (.not. associated(this%sigma_r)) return
    if (.not. allocated(this%Ml)) return

    select type(pGcache => this%G_r)
    type is (TMatrixCacheMem)
      if (.not. pGcache%isInitialized) return
      nE = size(pGcache%MatArray, 1)
      nK = size(pGcache%MatArray, 2)
    class default
      return
    end select

    nbl = this%struct%num_PLs
    pl = this%pl_index
    orb = this%orb_index
    pl_start = this%struct%mat_PL_start(pl)
    pl_end = this%struct%mat_PL_end(pl)
    nrow = pl_end - pl_start + 1

    select type(pS => this%sigma_r)
    type is (TMatrixCacheMem)
      if (.not. pS%isInitialized) then
        call pS%init(nE, nK, nbl, 3, 1)
      end if
    end select

    ! Capture into empty defect ~ (1-fd)*(Nq+1); emission from occupied ~ fd*Nq
    fac_em  = (this%Nq + 1.0_dp) * (1.0_dp - this%fd)
    fac_abs = this%Nq * this%fd

    label%spin = 1
    label%row_block = pl
    label%col_block = pl
    labelp = label

    do iK = 1, nK
      do iE = 1, nE
        label%kpoint = iK
        label%energy_point = iE
        acc = (0.0_dp, 0.0_dp)

        do l = 0, this%l_max
          wgt = this%Ml(l)
          if (wgt <= 0.0_dp) cycle

          iEm = iE - l * this%iEshift
          if (iEm >= 1 .and. iEm <= nE) then
            labelp%kpoint = iK
            labelp%energy_point = iEm
            if (this%G_r%is_cached(labelp)) then
              call this%G_r%retrieve_pointer(pG, labelp)
              if (associated(pG)) then
                if (allocated(pG%val)) then
                  if (orb <= size(pG%val, 1)) then
                    acc = acc + wgt * fac_em * pG%val(orb, orb)
                  end if
                end if
              end if
            end if
          end if

          if (l == 0) cycle
          iEp = iE + l * this%iEshift
          if (iEp >= 1 .and. iEp <= nE) then
            labelp%kpoint = iK
            labelp%energy_point = iEp
            if (this%G_r%is_cached(labelp)) then
              call this%G_r%retrieve_pointer(pG, labelp)
              if (associated(pG)) then
                if (allocated(pG%val)) then
                  if (orb <= size(pG%val, 1)) then
                    acc = acc + wgt * fac_abs * pG%val(orb, orb)
                  end if
                end if
              end if
            end if
          end if
        end do

        if (acc == (0.0_dp, 0.0_dp)) cycle

        call create(sig, nrow, nrow)
        sig%val = (0.0_dp, 0.0_dp)
        sig%val(orb, orb) = this%V2 * acc
        call this%sigma_r%add(sig, label)
        call mat_destroy(sig)
      end do
    end do

  end subroutine compute_Sigma_r

  subroutine compute_Sigma_n(this, en_index, k_index, spin)
    class(TSrhMultiphonon) :: this
    integer, intent(in), optional :: en_index, k_index, spin

    integer :: nE, nK, nbl, pl, orb, nrow, iE, iK, pl_start, pl_end, l, iEp, iEm
    type(TMatLabel) :: label, labelp
    type(z_dns), pointer :: pG
    type(z_dns) :: sig
    complex(dp) :: acc
    real(dp) :: fac_em, fac_abs, wgt

    if (.not. site_ok(this)) return
    if (this%V2 <= 0.0_dp) return
    if (.not. associated(this%G_n)) return
    if (.not. associated(this%sigma_n)) return
    if (.not. allocated(this%Ml)) return

    select type(pGcache => this%G_n)
    type is (TMatrixCacheMem)
      if (.not. pGcache%isInitialized) return
      nE = size(pGcache%MatArray, 1)
      nK = size(pGcache%MatArray, 2)
    class default
      return
    end select

    nbl = this%struct%num_PLs
    pl = this%pl_index
    orb = this%orb_index
    pl_start = this%struct%mat_PL_start(pl)
    pl_end = this%struct%mat_PL_end(pl)
    nrow = pl_end - pl_start + 1

    select type(pS => this%sigma_n)
    type is (TMatrixCacheMem)
      if (.not. pS%isInitialized) then
        call pS%init(nE, nK, nbl, 3, 1)
      end if
    end select

    ! Same fd weighting as Sigma_r (from previous SCBA iter / init fd=0.5)
    fac_em  = (this%Nq + 1.0_dp) * (1.0_dp - this%fd)
    fac_abs = this%Nq * this%fd

    label%spin = 1
    label%row_block = pl
    label%col_block = pl
    labelp = label

    do iK = 1, nK
      do iE = 1, nE
        label%kpoint = iK
        label%energy_point = iE
        acc = (0.0_dp, 0.0_dp)

        do l = 0, this%l_max
          wgt = this%Ml(l)
          if (wgt <= 0.0_dp) cycle

          iEm = iE - l * this%iEshift
          if (iEm >= 1 .and. iEm <= nE) then
            labelp%kpoint = iK
            labelp%energy_point = iEm
            if (this%G_n%is_cached(labelp)) then
              call this%G_n%retrieve_pointer(pG, labelp)
              if (associated(pG)) then
                if (allocated(pG%val)) then
                  if (orb <= size(pG%val, 1)) then
                    acc = acc + wgt * fac_em * pG%val(orb, orb)
                  end if
                end if
              end if
            end if
          end if

          if (l == 0) cycle
          iEp = iE + l * this%iEshift
          if (iEp >= 1 .and. iEp <= nE) then
            labelp%kpoint = iK
            labelp%energy_point = iEp
            if (this%G_n%is_cached(labelp)) then
              call this%G_n%retrieve_pointer(pG, labelp)
              if (associated(pG)) then
                if (allocated(pG%val)) then
                  if (orb <= size(pG%val, 1)) then
                    acc = acc + wgt * fac_abs * pG%val(orb, orb)
                  end if
                end if
              end if
            end if
          end if
        end do

        if (acc == (0.0_dp, 0.0_dp)) cycle

        call create(sig, nrow, nrow)
        sig%val = (0.0_dp, 0.0_dp)
        sig%val(orb, orb) = this%V2 * acc
        call this%sigma_n%add(sig, label)
        call mat_destroy(sig)
      end do
    end do

    call this%update_occupation()

  end subroutine compute_Sigma_n

  !> Stage 7: continuum-based defect occupation and net rate diagnostic.
  !! Stage-1 / classical two-carrier form (Aeberhard bulk limit):
  !!   n_loc from integral of local G^n at (pl,orb)  [NEGF-internal scale]
  !!   n1,p1 from Ed0 relative to midgap; ni diagnostic default
  !!   p_loc diagnostic: equal-injection proxy (= n_loc) until VB continuum is mapped
  !!   fd = (Cn*n + Cp*p1) / (Cn*(n+n1) + Cp*(p+p1))
  !!   r_net = (n*p - ni^2) / (tau_p*(n+n1) + tau_n*(p+p1))
  !! Electron-only balance would force r_net=0 by construction; do not use that.
  subroutine update_occupation(this)
    class(TSrhMultiphonon) :: this

    integer :: nE, nK, iE, iK, pl, orb
    type(TMatLabel) :: label
    type(z_dns), pointer :: pG
    real(dp) :: dE, kT, ni, n_sum, pi
    real(dp) :: n, p, n1, p1, tau_n, tau_p, denom, num

    if (.not. site_ok(this)) return
    if (.not. associated(this%G_n)) return

    select type(pGcache => this%G_n)
    type is (TMatrixCacheMem)
      if (.not. pGcache%isInitialized) return
      nE = size(pGcache%MatArray, 1)
      nK = size(pGcache%MatArray, 2)
    class default
      return
    end select

    pl = this%pl_index
    orb = this%orb_index
    dE = this%dE
    if (dE <= 0.0_dp) dE = 1.0_dp
    pi = 3.14159265358979323846_dp

    n_sum = 0.0_dp
    label%spin = 1
    label%row_block = pl
    label%col_block = pl

    do iK = 1, nK
      do iE = 1, nE
        label%kpoint = iK
        label%energy_point = iE
        if (.not. this%G_n%is_cached(label)) cycle
        call this%G_n%retrieve_pointer(pG, label)
        if (.not. associated(pG)) cycle
        if (.not. allocated(pG%val)) cycle
        if (orb > size(pG%val, 1)) cycle
        n_sum = n_sum + max(-aimag(pG%val(orb, orb)), 0.0_dp) * dE / (2.0_dp * pi)
      end do
    end do

    this%n_loc = n_sum / real(max(nK, 1), dp)
    n = this%n_loc
    ! Until VB continuum is available: equal-injection diagnostic proxy
    p = this%n_loc

    kT = KB_EV * max(this%temperature, 1.0_dp)
    ! ni on NEGF-internal scale. Factor 0.1 → mild injection so r_net is not
    ! identically zero (equal n=p=ni would force np=ni^2).
    ni = max(0.1_dp * this%n_loc, 1.0e-30_dp)
    n1 = ni * exp( this%Ed0 / kT)
    p1 = ni * exp(-this%Ed0 / kT)
    this%n1 = n1

    if (this%Cn * this%rho_d > 0.0_dp) then
      tau_n = 1.0_dp / (this%Cn * this%rho_d)
    else
      tau_n = 1.0e30_dp
    end if
    if (this%Cp * this%rho_d > 0.0_dp) then
      tau_p = 1.0_dp / (this%Cp * this%rho_d)
    else
      tau_p = 1.0e30_dp
    end if

    denom = this%Cn * (n + n1) + this%Cp * (p + p1)
    if (denom > 0.0_dp) then
      this%fd = (this%Cn * n + this%Cp * p1) / denom
    else
      this%fd = 0.5_dp
    end if

    denom = tau_p * (n + n1) + tau_n * (p + p1)
    num = n * p - ni * ni
    if (denom > 0.0_dp) then
      this%r_net = num / denom
    else
      this%r_net = 0.0_dp
    end if

    write(*,'(A,I3,A,ES12.4,A,F8.5,A,ES12.4)') &
        & ' SRH diagnostic (scba_iter=', this%scba_iter, &
        & '): n_loc=', this%n_loc, ' fd=', this%fd, ' r_net=', this%r_net

  end subroutine update_occupation

  subroutine destroy_Sigma_r(this)
    class(TSrhMultiphonon) :: this
    if (associated(this%sigma_r)) then
      call this%sigma_r%destroy()
    end if
  end subroutine destroy_Sigma_r

  subroutine destroy_Sigma_n(this)
    class(TSrhMultiphonon) :: this
    if (associated(this%sigma_n)) then
      call this%sigma_n%destroy()
    end if
  end subroutine destroy_Sigma_n

  subroutine destroy(this)
    class(TSrhMultiphonon) :: this
    call this%destroy_Sigma_r()
    call this%destroy_Sigma_n()
    if (associated(this%sigma_r)) then
      deallocate(this%sigma_r)
      nullify(this%sigma_r)
    end if
    if (associated(this%sigma_n)) then
      deallocate(this%sigma_n)
      nullify(this%sigma_n)
    end if
    if (allocated(this%Ml)) deallocate(this%Ml)
  end subroutine destroy

end module srh_multiphonon
