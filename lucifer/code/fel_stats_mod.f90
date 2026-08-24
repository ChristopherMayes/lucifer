!+
! Module fel_stats_mod
!
! The tracker's production statistics file <out_root>.stats.h5 (manual sec:stats).
!
! Fixed Bmad units everywhere (m, rad, eV, s, C, J, W), with units attributes written
! as documentation only, never load-bearing. Per-record datasets are (nrec, nslice)
! arrays as read by h5py (Genesis4-style per-slice export for visualization), with beam
! datasets named EXACTLY as bunch_params_struct components, sufficient to construct one
! from any (record, slice): centroid (nrec, nslice, 6), sigma (nrec, nslice, 36: the
! 6x6 flattened, symmetric, so the flattening order cannot mislead), charge_live,
! charge_tot, n_particle_live, n_particle_tot, t, sigma_t, plus z(nrec) and s == z.
! The field side stores wavefront_params_struct components per slice: centroid
! (nrec, nslice, 4), sigma (nrec, nslice, 16), energy, power, on_axis_intensity,
! emit_x, emit_y, angle_moments_valid (0/1). Theta moments cost FFTs, so they fill
! at element ends only (the twiss_valid pattern). Pulse-level values are POOLED
! downstream (scripts), never stored: the file stays raw.
!
! At ELEMENT ENDS the fully evaluated Bmad bunch_params_struct lands under
! element_end/: per whole-window bunch AND per slice, via Bmad's own calc_bunch_params
! (the Tao end-of-element pattern): twiss groups x/y/z/a/b/c each carrying beta,
! alpha, gamma, emit, norm_emit, sigma, sigma_p, eta, etap; plus centroid, sigma,
! charge_live, n_particle_live, ix_ele, s, twiss_valid. scripts/
! bunch_params_from_stats.py reconstructs a bunch_params dict from any (record, slice)
! of the per-record sufficient statistics. The harness checks that at element ends it
! reproduces these stored values.
!
! Beam moments are computed two-pass (mean first, then centered second moments, the
! manual sec:numerics variance rule: the one-pass form loses the entire sigma to
! cancellation), weighted by macroparticle charge, parallel over slices with
! fixed-order results (each slice's sums are its own).
!
! t and sigma_t derive from the stored chart: t = z_now/c - <z>/(beta0 c) and
! sigma_t = sigma_z/(beta0 c) with beta0 the reference beta. Exact for the reference
! particle, and documented rather than hidden inside a convention.
!-

module fel_stats_mod

use bmad
use beam_utils, only: calc_bunch_params, calc_emittances_and_twiss_from_sigma_matrix
use fel_beam_mod
use wavefront_mod
use fel_track_mod, only: fel_slip_struct, fel_field_struct, fel_field_index, fel_field_diag, fel_slice_to_bunch
use hdf5_interface

implicit none

!+
! Structure fel_stats_struct
!
! The accumulated statistics of one run: the per-record per-slice beam and field
! arrays, and the element-end evaluated bunch_params rows. Sized once by
! fel_stats_init from exact counts precomputed on the lattice walk.
!-

type fel_stats_struct
  real(rp) :: p0c = 0                 ! Reference momentum [eV] (for norm_emit reconstruction).
  integer :: nslice = 0
  integer :: nrec = 0, irec = 0       ! Capacity / fill of per-record arrays.
  integer :: nend = 0, iend = 0       ! Capacity / fill of element-end arrays.
  ! Per record, per slice. Beam side, bunch_params_struct names.
  real(rp), allocatable :: z(:)                   ! (nrec) [m]
  real(rp), allocatable :: b_centroid(:,:,:)      ! (6, nslice, nrec)
  real(rp), allocatable :: b_sigma(:,:,:)         ! (36, nslice, nrec)
  real(rp), allocatable :: charge_live(:,:)       ! (nslice, nrec) [C]
  real(rp), allocatable :: t(:,:), sigma_t(:,:)   ! (nslice, nrec) [s]
  real(rp), allocatable :: bunching(:,:)          ! (nslice, nrec) |b| at the fundamental
  real(rp), allocatable :: bunching_phase(:,:)    ! (nslice, nrec) [rad]
  integer, allocatable :: n_particle_live(:,:)    ! (nslice, nrec)
  ! Field side, wavefront_params_struct names. With ONE live polarization these are
  ! the whole story. With two, they carry the X component, the f2_* arrays carry Y
  ! (written as a field/y/ group), and the power/energy/on_axis datasets are written
  ! as TOTALS. Single-polarization files are unchanged.
  real(rp), allocatable :: f_centroid(:,:,:)      ! (4, nslice, nrec)
  real(rp), allocatable :: f_sigma(:,:,:)         ! (16, nslice, nrec)
  real(rp), allocatable :: f_energy(:,:), f_power(:,:), f_on_axis(:,:)
  real(rp), allocatable :: f_emit_x(:,:), f_emit_y(:,:)
  integer, allocatable :: f_angles_valid(:,:)     ! 0/1
  real(rp), allocatable :: f2_centroid(:,:,:), f2_sigma(:,:,:)
  real(rp), allocatable :: f2_energy(:,:), f2_power(:,:), f2_on_axis(:,:)
  real(rp), allocatable :: f2_emit_x(:,:), f2_emit_y(:,:)
  ! Harmonic fields (field-set entries 2+): full wavefront_params per harmonic,
  ! written as field/harm<h>/ groups. The fundamental's datasets above are NEVER
  ! summed with these: harmonics are distinct wavelengths (a detector separates
  ! colors), unlike the two polarization components of one wave.
  integer, allocatable :: fh_harm(:)                ! (n_extra) harmonic numbers.
  real(rp), allocatable :: fh_centroid(:,:,:,:)     ! (4, nslice, nrec, n_extra)
  real(rp), allocatable :: fh_sigma(:,:,:,:)        ! (16, nslice, nrec, n_extra)
  real(rp), allocatable :: fh_energy(:,:,:), fh_power(:,:,:), fh_on_axis(:,:,:)
  real(rp), allocatable :: fh_emit_x(:,:,:), fh_emit_y(:,:,:)
  integer, allocatable :: fh_angles_valid(:,:,:)
  ! Element ends: the evaluated bunch_params_struct, whole bunch and per slice.
  integer, allocatable :: e_ix_ele(:)             ! (nend)
  real(rp), allocatable :: e_s(:)                 ! (nend) [m]
  real(rp), allocatable :: e_f_power(:,:)         ! (nslice, nend) [W] radiation power, TOTAL over
                                                  !   polarizations, as the per-record dataset is.
  real(rp), allocatable :: e_f_energy(:,:)        ! (nslice, nend) [J] total.
  real(rp), allocatable :: e_f_on_axis(:,:)       ! (nslice, nend) total.
  real(rp), allocatable :: e_bunching(:,:)        ! (nslice, nend) |b| at the fundamental.
  real(rp), allocatable :: e_bunching_phase(:,:)  ! (nslice, nend) [rad]
  real(rp), allocatable :: e_bunch(:,:)           ! (n_bp, nend)
  real(rp), allocatable :: e_slice(:,:,:)         ! (n_bp, nslice, nend)
end type

! The flattened bunch_params record: 6 centroid + 36 sigma + charge_live +
! n_particle_live + twiss_valid + 6 modes x 9 params = 99.
integer, parameter :: fel_stats_n_bp$ = 99
character(*), parameter :: fel_stats_modes$(6) = [character(1):: 'x', 'y', 'z', 'a', 'b', 'c']
character(9), parameter :: fel_stats_twiss$(9) = [character(9):: &
        'beta', 'alpha', 'gamma', 'emit', 'norm_emit', 'sigma', 'sigma_p', 'eta', 'etap']

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_stats_init (stats, nrec, nend, nslice, p0c, two_pol, harm_extra)
!
! Routine to size the accumulation arrays. nrec and nend are exact counts precomputed
! from the lattice walk (records = 1 + sum of undulator steps + one per interlude).
!
! Input:
!   nrec           -- integer: Number of per-slice records to allocate.
!   nend           -- integer: Number of element-end rows to allocate.
!   nslice         -- integer: Number of beam slices.
!   p0c            -- real(rp): Reference momentum [eV] (for norm_emit reconstruction).
!   two_pol        -- logical: If True, allocate the second-polarization (f2_*) arrays.
!   harm_extra(:)  -- integer: Harmonic numbers of field-set entries 2+ (may be empty).
!
! Output:
!   stats          -- fel_stats_struct: Allocated accumulator with zero fill counts.
!-

subroutine fel_stats_init (stats, nrec, nend, nslice, p0c, two_pol, harm_extra)

type (fel_stats_struct) stats
real(rp) p0c
integer nrec, nend, nslice
integer harm_extra(:)
logical two_pol

!

stats%p0c = p0c
stats%nslice = nslice
stats%nrec = nrec;  stats%irec = 0
stats%nend = nend;  stats%iend = 0

allocate (stats%z(nrec))
allocate (stats%b_centroid(6, nslice, nrec), stats%b_sigma(36, nslice, nrec))
allocate (stats%charge_live(nslice, nrec), stats%t(nslice, nrec), stats%sigma_t(nslice, nrec))
allocate (stats%bunching(nslice, nrec), stats%bunching_phase(nslice, nrec))
allocate (stats%n_particle_live(nslice, nrec))
allocate (stats%f_centroid(4, nslice, nrec), stats%f_sigma(16, nslice, nrec))
allocate (stats%f_energy(nslice, nrec), stats%f_power(nslice, nrec), stats%f_on_axis(nslice, nrec))
allocate (stats%f_emit_x(nslice, nrec), stats%f_emit_y(nslice, nrec))
allocate (stats%f_angles_valid(nslice, nrec))
if (two_pol) then
  allocate (stats%f2_centroid(4, nslice, nrec), stats%f2_sigma(16, nslice, nrec))
  allocate (stats%f2_energy(nslice, nrec), stats%f2_power(nslice, nrec), stats%f2_on_axis(nslice, nrec))
  allocate (stats%f2_emit_x(nslice, nrec), stats%f2_emit_y(nslice, nrec))
endif
if (size(harm_extra) > 0) then
  allocate (stats%fh_harm(size(harm_extra)))
  stats%fh_harm = harm_extra
  allocate (stats%fh_centroid(4, nslice, nrec, size(harm_extra)), &
            stats%fh_sigma(16, nslice, nrec, size(harm_extra)))
  allocate (stats%fh_energy(nslice, nrec, size(harm_extra)), &
            stats%fh_power(nslice, nrec, size(harm_extra)), &
            stats%fh_on_axis(nslice, nrec, size(harm_extra)))
  allocate (stats%fh_emit_x(nslice, nrec, size(harm_extra)), &
            stats%fh_emit_y(nslice, nrec, size(harm_extra)))
  allocate (stats%fh_angles_valid(nslice, nrec, size(harm_extra)))
endif
allocate (stats%e_ix_ele(nend), stats%e_s(nend))
allocate (stats%e_bunch(fel_stats_n_bp$, nend), stats%e_slice(fel_stats_n_bp$, nslice, nend))
allocate (stats%e_f_power(nslice, nend), stats%e_f_energy(nslice, nend), stats%e_f_on_axis(nslice, nend))
allocate (stats%e_bunching(nslice, nend), stats%e_bunching_phase(nslice, nend))
stats%e_f_power = 0;  stats%e_f_energy = 0;  stats%e_f_on_axis = 0
stats%e_bunching = 0; stats%e_bunching_phase = 0

end subroutine fel_stats_init

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_stats_record (stats, beam, ff, z_now, with_angles, bdiag_arr, fpow, fonax, err_flag)
!
! Routine to take one per-record row: beam sufficient statistics and wavefront params
! for every slice. with_angles fills the field theta moments (element ends). The same
! loop also evaluates the diag instrument for every slice, so the diag writer then
! only prints.
!
! Input:
!   stats         -- fel_stats_struct: Accumulator.
!   beam          -- fel_beam_struct: The sliced beam.
!   ff(:)         -- fel_field_struct: The field set. Entry 1 is the fundamental.
!   z_now         -- real(rp): Position of the record [m].
!   with_angles   -- logical: If True fill the field theta moments (they cost FFTs,
!                      so element ends only).
!
! Output:
!   stats         -- fel_stats_struct: Row stats%irec filled.
!   bdiag_arr(:)  -- fel_slice_diag_struct: The diag instrument, one per slice.
!   fpow(:)       -- real(rp): fel_field_diag's power per slice [W].
!   fonax(:)      -- real(rp): fel_field_diag's on-axis intensity per slice.
!   err_flag      -- logical: Set True on a wavefront_params error or record overflow.
!-

subroutine fel_stats_record (stats, beam, ff, z_now, with_angles, bdiag_arr, fpow, fonax, err_flag)

type (fel_stats_struct) stats
type (fel_beam_struct), target :: beam
type (fel_field_struct), target :: ff(:)
type (wavefront_struct), pointer :: wf
type (fel_slice_struct), pointer :: sl
type (wavefront_params_struct) pms
integer io
type (fel_slice_diag_struct) bdiag_arr(:)
real(rp) fpow(:), fonax(:)
real(rp) z_now
logical with_angles, err_flag

real(rp) w, wsum, mean(6), cen(6), sig(6,6), v(6), beta0, ks
integer ir, is, ip, i, j, nslice
logical err, any_err
character(*), parameter :: r_name = 'fel_stats_record'

!

err_flag = .true.
wf => ff(1)%wf
nslice = size(beam%slice)
if (stats%irec >= stats%nrec) then
  call out_io (s_error$, r_name, 'MORE RECORDS THAN THE PRECOMPUTED COUNT. INTERNAL BUG.')
  return
endif
stats%irec = stats%irec + 1
ir = stats%irec
stats%z(ir) = z_now
beta0 = fel_p0_mc(beam) / sqrt(fel_p0_mc(beam)**2 + 1)
ks = twopi / wf%wavelength

any_err = .false.

! This loop ALSO evaluates the diag instrument (fel_slice_diag, fel_field_diag) for
! every slice: the diag writer then only prints. Each slice's arithmetic is the
! identical serial code, so diag.txt is bit-for-bit what it always was. What changed
! is that the formerly SERIAL per-record diag sweeps now ride this parallel loop.

!$OMP parallel do private(sl, w, wsum, mean, cen, sig, v, ip, i, j, io, pms, err) &
!$OMP&   reduction(.or.: any_err)
do is = 1, nslice
  sl => beam%slice(is)
  call fel_field_diag (wf, fel_field_index(ff(1)%slip, is, nslice), fpow(is), fonax(is))
  call fel_slice_diag (beam, sl, ks, bdiag_arr(is))

  ! Two-pass weighted moments: mean first, then centered second moments (manual
  ! sec:numerics).

  wsum = 0;  mean = 0
  do ip = 1, sl%n
    w = sl%weight(ip)
    wsum = wsum + w
    mean = mean + w * [sl%x(ip), sl%px(ip), sl%y(ip), sl%py(ip), sl%z(ip), sl%pz(ip)]
  enddo
  if (wsum > 0) mean = mean / wsum

  sig = 0
  do ip = 1, sl%n
    v = [sl%x(ip), sl%px(ip), sl%y(ip), sl%py(ip), sl%z(ip), sl%pz(ip)] - mean
    w = sl%weight(ip)
    do j = 1, 6
      do i = 1, j
        sig(i,j) = sig(i,j) + w * v(i) * v(j)
      enddo
    enddo
  enddo
  if (wsum > 0) sig = sig / wsum
  do j = 1, 6
    do i = j+1, 6
      sig(i,j) = sig(j,i)
    enddo
  enddo

  cen = mean
  stats%b_centroid(:, is, ir) = cen
  stats%b_sigma(:, is, ir) = reshape(sig, [36])
  stats%charge_live(is, ir) = wsum
  stats%n_particle_live(is, ir) = sl%n
  stats%t(is, ir) = z_now / c_light - cen(5) / (beta0 * c_light)
  stats%sigma_t(is, ir) = sqrt(max(0.0_rp, sig(5,5))) / (beta0 * c_light)

  ! Bunching from the diag instrument just evaluated: |b| at the fundamental. The
  ! phase carries the common phi0, Genesis's own convention.

  stats%bunching(is, ir) = bdiag_arr(is)%bunching
  stats%bunching_phase(is, ir) = bdiag_arr(is)%bunching_phase

  ! The field slice this beam slice couples to, unrotated exactly as the dumps are.

  call wavefront_params_of_plane (wf%Ex(:,:,fel_field_index(ff(1)%slip, is, nslice)), wf%dx, &
                                  wf%wavelength, beam%slice_spacing, pms, with_angles, err)
  any_err = any_err .or. err
  pms%s = z_now
  stats%f_centroid(:, is, ir) = pms%centroid
  stats%f_sigma(:, is, ir) = reshape(pms%sigma, [16])
  stats%f_energy(is, ir) = pms%energy
  stats%f_power(is, ir) = pms%power
  stats%f_on_axis(is, ir) = pms%on_axis_intensity
  stats%f_emit_x(is, ir) = pms%emit_x
  stats%f_emit_y(is, ir) = pms%emit_y
  stats%f_angles_valid(is, ir) = merge(1, 0, pms%angle_moments_valid)

  if (allocated(wf%Ey)) then
    call wavefront_params_of_plane (wf%Ey(:,:,fel_field_index(ff(1)%slip, is, nslice)), wf%dx, &
                                    wf%wavelength, beam%slice_spacing, pms, with_angles, err)
    any_err = any_err .or. err
    stats%f2_centroid(:, is, ir) = pms%centroid
    stats%f2_sigma(:, is, ir) = reshape(pms%sigma, [16])
    stats%f2_energy(is, ir) = pms%energy
    stats%f2_power(is, ir) = pms%power
    stats%f2_on_axis(is, ir) = pms%on_axis_intensity
    stats%f2_emit_x(is, ir) = pms%emit_x
    stats%f2_emit_y(is, ir) = pms%emit_y
  endif

  ! Harmonic fields: full wavefront_params per extra field, at ITS wavelength and
  ! record rotation.

  do io = 2, size(ff)
    call wavefront_params_of_plane ( &
            ff(io)%wf%Ex(:,:,fel_field_index(ff(io)%slip, is, size(ff(io)%wf%Ex,3))), &
            ff(io)%wf%dx, ff(io)%wf%wavelength, beam%slice_spacing, pms, with_angles, err)
    any_err = any_err .or. err
    stats%fh_centroid(:, is, ir, io-1) = pms%centroid
    stats%fh_sigma(:, is, ir, io-1) = reshape(pms%sigma, [16])
    stats%fh_energy(is, ir, io-1) = pms%energy
    stats%fh_power(is, ir, io-1) = pms%power
    stats%fh_on_axis(is, ir, io-1) = pms%on_axis_intensity
    stats%fh_emit_x(is, ir, io-1) = pms%emit_x
    stats%fh_emit_y(is, ir, io-1) = pms%emit_y
    stats%fh_angles_valid(is, ir, io-1) = merge(1, 0, pms%angle_moments_valid)
  enddo
enddo
!$OMP end parallel do

err_flag = any_err

end subroutine fel_stats_record

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_stats_element_end (stats, beam, ele, z_now, err_flag)
!
! Routine to take one element-end row: Bmad's calc_bunch_params for the whole window
! (all slices as one bunch in global window coordinates) and per slice (the Tao
! end-of-element pattern), using Bmad's own statistics machinery.
!
! The row is SELF-SUFFICIENT: beam moments and Twiss for the whole window and per
! slice, plus the radiation power, energy, on-axis intensity and bunching per slice.
! That matters when comb_ds_save < 0 keeps no per-record rows at all (Bmad's comb
! semantics, kept verbatim). The element-end row is then the only record of the run,
! and it carries the field as well as the beam.
!
! Input:
!   stats     -- fel_stats_struct: Accumulator. The current record supplies the slice moments.
!   beam      -- fel_beam_struct: The sliced beam.
!   ff(:)     -- fel_field_struct: The field set, for the field values when there is no record.
!   ele       -- ele_struct: The element just ended.
!   z_now     -- real(rp): Position of the element end [m].
!
! Output:
!   stats     -- fel_stats_struct: Row stats%iend filled.
!   err_flag  -- logical: Set True on a conversion error or row overflow. False otherwise.
!-

subroutine fel_stats_element_end (stats, beam, ff, ele, z_now, err_flag)

type (fel_stats_struct) stats
type (fel_beam_struct), target :: beam
type (fel_field_struct), target :: ff(:)
type (wavefront_struct), pointer :: wf
type (fel_slice_struct), pointer :: sl
type (wavefront_params_struct) pms
type (fel_slice_diag_struct) bdg
type (ele_struct) ele
type (bunch_params_struct) bp
real(rp), allocatable :: cen_s(:,:), sig_s(:,:,:), w_s(:), dz_s(:), shear_s(:)
integer, allocatable :: n_s(:)
real(rp) z_now, p0_mc, p_mc, w, wsum, mean(6), sig(6,6), v(6), ks
integer ie, is, ip, i, j, nslice
logical err_flag, error, err, any_err
character(*), parameter :: r_name = 'fel_stats_element_end'

!

err_flag = .true.
if (stats%iend >= stats%nend) then
  call out_io (s_error$, r_name, 'MORE ELEMENT ENDS THAN THE PRECOMPUTED COUNT. INTERNAL BUG.')
  return
endif
stats%iend = stats%iend + 1
ie = stats%iend
stats%e_ix_ele(ie) = ele%ix_ele
stats%e_s(ie) = z_now

nslice = size(beam%slice)
wf => ff(1)%wf
allocate (cen_s(6,nslice), sig_s(6,6,nslice), w_s(nslice), dz_s(nslice), shear_s(nslice), n_s(nslice))

! The field and bunching side of the row. With a current record they are already
! evaluated there, so copy (the datasets then agree exactly). With no records evaluate
! them here, the same routines the record sweep uses, angle moments not needed.

if (stats%irec > 0) then
  stats%e_f_power(:,ie) = stats%f_power(:, stats%irec)
  stats%e_f_energy(:,ie) = stats%f_energy(:, stats%irec)
  stats%e_f_on_axis(:,ie) = stats%f_on_axis(:, stats%irec)
  if (allocated(stats%f2_energy)) then
    stats%e_f_power(:,ie) = stats%e_f_power(:,ie) + stats%f2_power(:, stats%irec)
    stats%e_f_energy(:,ie) = stats%e_f_energy(:,ie) + stats%f2_energy(:, stats%irec)
    stats%e_f_on_axis(:,ie) = stats%e_f_on_axis(:,ie) + stats%f2_on_axis(:, stats%irec)
  endif
  stats%e_bunching(:,ie) = stats%bunching(:, stats%irec)
  stats%e_bunching_phase(:,ie) = stats%bunching_phase(:, stats%irec)

else
  any_err = .false.
  ks = twopi / wf%wavelength
  !$OMP parallel do private(sl, pms, bdg, err) reduction(.or.: any_err)
  do is = 1, nslice
    sl => beam%slice(is)
    call wavefront_params_of_plane (wf%Ex(:,:,fel_field_index(ff(1)%slip, is, nslice)), wf%dx, &
                                    wf%wavelength, beam%slice_spacing, pms, .false., err)
    any_err = any_err .or. err
    stats%e_f_power(is,ie) = pms%power
    stats%e_f_energy(is,ie) = pms%energy
    stats%e_f_on_axis(is,ie) = pms%on_axis_intensity
    if (allocated(wf%Ey)) then
      call wavefront_params_of_plane (wf%Ey(:,:,fel_field_index(ff(1)%slip, is, nslice)), wf%dx, &
                                      wf%wavelength, beam%slice_spacing, pms, .false., err)
      any_err = any_err .or. err
      stats%e_f_power(is,ie) = stats%e_f_power(is,ie) + pms%power
      stats%e_f_energy(is,ie) = stats%e_f_energy(is,ie) + pms%energy
      stats%e_f_on_axis(is,ie) = stats%e_f_on_axis(is,ie) + pms%on_axis_intensity
    endif
    call fel_slice_diag (beam, sl, ks, bdg)
    stats%e_bunching(is,ie) = bdg%bunching
    stats%e_bunching_phase(is,ie) = bdg%bunching_phase
  enddo
  !$OMP end parallel do
  if (any_err) return
endif

! ONE source of per-slice moments for both the per-slice rows and the whole-window
! row. With per-record rows (the comb >= 0) the current record already holds them:
! an element end always coincides with its last record. With NO per-record rows (the
! comb's "< 0 => No comb calculated") take them here, in one parallel sweep whose
! two-pass weighted arithmetic is the per-record sweep's own.

if (stats%irec > 0) then
  do is = 1, nslice
    cen_s(:,is) = stats%b_centroid(:, is, stats%irec)
    sig_s(:,:,is) = reshape(stats%b_sigma(:, is, stats%irec), [6, 6])
    w_s(is) = stats%charge_live(is, stats%irec)
    n_s(is) = stats%n_particle_live(is, stats%irec)
  enddo

else
  !$OMP parallel do private(sl, ip, i, j, w, wsum, mean, sig, v)
  do is = 1, nslice
    sl => beam%slice(is)
    wsum = 0;  mean = 0
    do ip = 1, sl%n
      w = sl%weight(ip)
      wsum = wsum + w
      mean = mean + w * [sl%x(ip), sl%px(ip), sl%y(ip), sl%py(ip), sl%z(ip), sl%pz(ip)]
    enddo
    if (wsum > 0) mean = mean / wsum

    sig = 0
    do ip = 1, sl%n
      v = [sl%x(ip), sl%px(ip), sl%y(ip), sl%py(ip), sl%z(ip), sl%pz(ip)] - mean
      w = sl%weight(ip)
      do j = 1, 6
        do i = 1, j
          sig(i,j) = sig(i,j) + w * v(i) * v(j)
        enddo
      enddo
    enddo
    if (wsum > 0) sig = sig / wsum
    do j = 1, 6
      do i = j+1, 6
        sig(i,j) = sig(j,i)
      enddo
    enddo

    cen_s(:,is) = mean;  sig_s(:,:,is) = sig
    w_s(is) = wsum;      n_s(is) = sl%n
  enddo
  !$OMP end parallel do
endif

! Per slice: the twiss evaluation is Bmad's own (calc_emittances_and_twiss_from_
! sigma_matrix, the identical code path calc_bunch_params ends in), fed from the
! moments above instead of re-summing every particle through a bunch conversion.
! Twiss only where there are enough live particles: a time window's near-empty edge
! slices are degenerate BY CONSTRUCTION (zero charge, collapsed sigma modes).

do is = 1, nslice
  bp = bunch_params_struct()
  bp%centroid%vec = cen_s(:,is)
  bp%centroid%p0c = stats%p0c
  bp%centroid%species = electron$
  bp%sigma = sig_s(:,:,is)
  bp%charge_live = w_s(is)
  bp%n_particle_live = n_s(is)
  bp%n_particle_tot = bp%n_particle_live
  if (bp%n_particle_live >= 6) then
    call calc_emittances_and_twiss_from_sigma_matrix (bp%sigma, bp, error, .false.)
  endif
  call pack_bp (bp, stats%e_slice(:, is, ie))
enddo

! The whole window, from the same per-slice moments, no particle visit at all. Each
! slice's stored moments live in its LOCAL z chart and enter the pool moved to the
! global window chart by the migration invariant fel_concat_slices uses,
! z_global = z_local + beta*(is-1)*slice_spacing.
!
! That map is NOT a constant offset: fel_concat_slices evaluates beta per particle, so
! within a slice z_global depends on pz. Linearizing beta about the slice's own mean pz
! makes the map an exact SHEAR of (z, pz),
!   z_global = z_local + L*beta_bar + (L*dbeta/dpz)*(pz - pz_bar),   L = (is-1)*spacing,
! whose effect on the slice covariance is S -> J S J^T with the single off-diagonal
! J(5,6) = L*dbeta/dpz. Dropping that shear leaves var(z) and cov(z,pz) right (they are
! dominated by the window-scale spread) but biases every cov(v,z) cross term, which are
! small by near-cancellation: measured 4e-6 relative on those terms, against 2e-13 for
! the rest, and it flipped the whole-window normal-mode decomposition from valid to
! invalid. The residual after the shear is second order, O(d2beta/dpz2 * sigma_pz^2).

p0_mc = fel_p0_mc(beam)
do is = 1, nslice
  p_mc = p0_mc * (1 + cen_s(6,is))
  dz_s(is) = fel_beta_of(p0_mc, cen_s(6,is)) * (is-1) * beam%slice_spacing
  shear_s(is) = (is-1) * beam%slice_spacing * p0_mc / sqrt(p_mc**2 + 1)**3
enddo

call fel_pool_bunch_params (cen_s, sig_s, w_s, n_s, dz_s, shear_s, stats%p0c, bp)
if (bp%n_particle_live >= 6) then
  call calc_emittances_and_twiss_from_sigma_matrix (bp%sigma, bp, error, .false.)
endif
call pack_bp (bp, stats%e_bunch(:, ie))

err_flag = .false.

!------------------------------------------------------------------------------
contains

!+
! Subroutine pack_bp (bp, row)
!
! Routine to flatten one bunch_params_struct into a stats row: centroid, sigma,
! charge_live, n_particle_live, twiss_valid, then the nine twiss parameters of each
! of the six modes.
!-

subroutine pack_bp (bp, row)

type (bunch_params_struct) bp
type (twiss_struct) tw
real(rp) row(:)
integer im, k

!

row(1:6) = bp%centroid%vec
row(7:42) = reshape(bp%sigma, [36])
row(43) = bp%charge_live
row(44) = bp%n_particle_live
row(45) = merge(1.0_rp, 0.0_rp, bp%twiss_valid)
k = 45
do im = 1, 6
  select case (fel_stats_modes$(im))
  case ('x'); tw = bp%x
  case ('y'); tw = bp%y
  case ('z'); tw = bp%z
  case ('a'); tw = bp%a
  case ('b'); tw = bp%b
  case ('c'); tw = bp%c
  end select
  row(k+1:k+9) = [tw%beta, tw%alpha, tw%gamma, tw%emit, tw%norm_emit, tw%sigma, tw%sigma_p, tw%eta, tw%etap]
  k = k + 9
enddo

end subroutine pack_bp

end subroutine fel_stats_element_end

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_pool_bunch_params (cen_s, sig_s, w_s, n_s, dz_s, p0c, bp)
!
! Routine to assemble one bunch_params_struct for the whole time window from the
! per-slice moments, without visiting particles. This is the beam-side twin of
! pool_wavefront (wavefront_mod), and it is the pooled-covariance identity, exact for
! any set of groups:
!
!   m = Sum_s w_s m_s / Sum_s w_s
!   S = Sum_s w_s [S_s + (m_s - m)(m_s - m)^T] / Sum_s w_s
!
! The between-group term (m_s - m)(m_s - m)^T is what makes it exact rather than an
! average of covariances: without it the window's sigma(5,5) would collapse from the
! window length squared to a slice length squared.
!
! Only the components pack_bp writes are filled (centroid, sigma, charge_live,
! n_particle_live and, by the caller, the twiss groups from the assembled sigma).
! rel_max/rel_min are extrema, beyond this identity's reach, and are not written to
! the stats file. They are deliberately left at their defaults.
!
! Input:
!   cen_s(:,:)    -- real(rp): (6, nslice) per-slice centroid, in each slice's local chart.
!   sig_s(:,:,:)  -- real(rp): (6, 6, nslice) per-slice covariance.
!   w_s(:)        -- real(rp): (nslice) per-slice live charge, the pooling weight [C].
!   n_s(:)        -- integer: (nslice) per-slice live particle count.
!   dz_s(:)       -- real(rp): (nslice) local-to-global-chart z offset of each slice [m].
!   shear_s(:)    -- real(rp): (nslice) dz_global/dpz of that map, J(5,6) of the shear.
!   p0c           -- real(rp): Reference momentum [eV], for the norm_emit scale.
!
! Output:
!   bp            -- bunch_params_struct: Whole-window centroid, sigma, charge and count.
!-

subroutine fel_pool_bunch_params (cen_s, sig_s, w_s, n_s, dz_s, shear_s, p0c, bp)

type (bunch_params_struct) bp
real(rp) cen_s(:,:), sig_s(:,:,:), w_s(:), dz_s(:), shear_s(:), p0c
integer n_s(:)
real(rp) wtot, m(6), c(6), d(6), s6(6,6), k
integer is, i, j, nslice

!

nslice = size(w_s)
wtot = sum(w_s)

bp = bunch_params_struct()
bp%centroid%p0c = p0c
bp%centroid%species = electron$
bp%n_particle_live = sum(n_s)
bp%n_particle_tot = bp%n_particle_live
bp%charge_live = wtot
bp%charge_tot = wtot
if (wtot <= 0) return

m = 0
do is = 1, nslice
  c = cen_s(:,is);  c(5) = c(5) + dz_s(is)
  m = m + w_s(is) * c
enddo
m = m / wtot

do is = 1, nslice
  c = cen_s(:,is);  c(5) = c(5) + dz_s(is)
  d = c - m

  ! S -> J S J^T with the one off-diagonal J(5,6) = shear_s: the z row and column pick
  ! up the pz coupling of the local-to-global map, in that order (the (5,5) update uses
  ! the ORIGINAL (5,6) and (6,6), so it goes first).

  s6 = sig_s(:,:,is)
  k = shear_s(is)
  if (k /= 0) then
    s6(5,5) = s6(5,5) + 2 * k * sig_s(5,6,is) + k**2 * sig_s(6,6,is)
    do j = 1, 6
      if (j == 5) cycle
      s6(5,j) = sig_s(5,j,is) + k * sig_s(6,j,is)
      s6(j,5) = s6(5,j)
    enddo
  endif

  do j = 1, 6
    do i = 1, 6
      bp%sigma(i,j) = bp%sigma(i,j) + w_s(is) * (s6(i,j) + d(i) * d(j))
    enddo
  enddo
enddo
bp%sigma = bp%sigma / wtot
bp%centroid%vec = m

end subroutine fel_pool_bunch_params

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_stats_write (stats, file_name, err_flag)
!
! Routine to write the stats file. Datasets appear to h5py with the record index first:
! Fortran (a, nslice, nrec) reads as (nrec, nslice, a). Units attributes are written as
! documentation. The units are FIXED and the attributes are never load-bearing.
!
! Input:
!   stats      -- fel_stats_struct: The filled accumulator.
!   file_name  -- character(*): File to write.
!
! Output:
!   err_flag   -- logical: Set True on a write error. False otherwise.
!-

subroutine fel_stats_write (stats, file_name, err_flag)

type (fel_stats_struct) stats
integer(hid_t) f_id, g_id, b_id, s_id
integer h5_err, ir, ihh
logical err_flag, err
character(*) file_name
character(12) hname
character(*), parameter :: r_name = 'fel_stats_write'

!

err_flag = .true.
ir = stats%irec

call hdf5_open_file (file_name, 'WRITE', f_id, err);  if (err) return

call hdf5_write_dataset_real (f_id, 'z', stats%z(1:ir), err);  if (err) return
call hdf5_write_dataset_real (f_id, 'p0c', [stats%p0c], err);  if (err) return

! One root attribute documents the fixed units. It is never load-bearing.

call hdf5_write_attribute_string (f_id, 'units_note', &
        'Fixed Bmad units: m, rad, eV, s, C, J, W. Bmad phase space (x, px/p0, y, py/p0, z, pz).', err)

call H5Gcreate_f (f_id, 'beam', g_id, h5_err)
call hdf5_write_dataset_real (g_id, 'centroid', stats%b_centroid(:,:,1:ir), err);  if (err) return
call hdf5_write_dataset_real (g_id, 'sigma', stats%b_sigma(:,:,1:ir), err);  if (err) return
call hdf5_write_dataset_real (g_id, 'charge_live', stats%charge_live(:,1:ir), err);  if (err) return
call hdf5_write_dataset_real (g_id, 'charge_tot', stats%charge_live(:,1:ir), err);  if (err) return
call hdf5_write_dataset_int (g_id, 'n_particle_live', stats%n_particle_live(:,1:ir), err);  if (err) return
call hdf5_write_dataset_int (g_id, 'n_particle_tot', stats%n_particle_live(:,1:ir), err);  if (err) return
call hdf5_write_dataset_real (g_id, 't', stats%t(:,1:ir), err);  if (err) return
call hdf5_write_dataset_real (g_id, 'sigma_t', stats%sigma_t(:,1:ir), err);  if (err) return
call hdf5_write_dataset_real (g_id, 'bunching', stats%bunching(:,1:ir), err);  if (err) return
call hdf5_write_dataset_real (g_id, 'bunching_phase', stats%bunching_phase(:,1:ir), err);  if (err) return
call hdf5_write_dataset_real (g_id, 's', stats%z(1:ir), err);  if (err) return
call H5Gclose_f (g_id, h5_err)

call H5Gcreate_f (f_id, 'field', g_id, h5_err)
call hdf5_write_dataset_real (g_id, 'centroid', stats%f_centroid(:,:,1:ir), err);  if (err) return
call hdf5_write_dataset_real (g_id, 'sigma', stats%f_sigma(:,:,1:ir), err);  if (err) return
if (allocated(stats%f2_energy)) then
  ! Two live polarizations: power/energy/intensity are TOTALS. The x-component params
  ! stay in this group's centroid/sigma/emit, the y component's under field/y/.
  call hdf5_write_dataset_real (g_id, 'energy', stats%f_energy(:,1:ir) + stats%f2_energy(:,1:ir), err);  if (err) return
  call hdf5_write_dataset_real (g_id, 'power', stats%f_power(:,1:ir) + stats%f2_power(:,1:ir), err);  if (err) return
  call hdf5_write_dataset_real (g_id, 'on_axis_intensity', &
                                stats%f_on_axis(:,1:ir) + stats%f2_on_axis(:,1:ir), err);  if (err) return
else
  call hdf5_write_dataset_real (g_id, 'energy', stats%f_energy(:,1:ir), err);  if (err) return
  call hdf5_write_dataset_real (g_id, 'power', stats%f_power(:,1:ir), err);  if (err) return
  call hdf5_write_dataset_real (g_id, 'on_axis_intensity', stats%f_on_axis(:,1:ir), err);  if (err) return
endif
call hdf5_write_dataset_real (g_id, 'emit_x', stats%f_emit_x(:,1:ir), err);  if (err) return
call hdf5_write_dataset_real (g_id, 'emit_y', stats%f_emit_y(:,1:ir), err);  if (err) return
call hdf5_write_dataset_int (g_id, 'angle_moments_valid', stats%f_angles_valid(:,1:ir), err);  if (err) return
if (allocated(stats%f2_energy)) then
  call H5Gcreate_f (g_id, 'y', b_id, h5_err)
  call hdf5_write_dataset_real (b_id, 'centroid', stats%f2_centroid(:,:,1:ir), err);  if (err) return
  call hdf5_write_dataset_real (b_id, 'sigma', stats%f2_sigma(:,:,1:ir), err);  if (err) return
  call hdf5_write_dataset_real (b_id, 'energy', stats%f2_energy(:,1:ir), err);  if (err) return
  call hdf5_write_dataset_real (b_id, 'power', stats%f2_power(:,1:ir), err);  if (err) return
  call hdf5_write_dataset_real (b_id, 'on_axis_intensity', stats%f2_on_axis(:,1:ir), err);  if (err) return
  call hdf5_write_dataset_real (b_id, 'emit_x', stats%f2_emit_x(:,1:ir), err);  if (err) return
  call hdf5_write_dataset_real (b_id, 'emit_y', stats%f2_emit_y(:,1:ir), err);  if (err) return
  call H5Gclose_f (b_id, h5_err)
endif

! Harmonic fields: field/harm<h>/ groups, one full wavefront_params set each. The
! fundamental's datasets above are untouched: harmonics are separate wavelengths.

if (allocated(stats%fh_energy)) then
  do ihh = 1, size(stats%fh_harm)
    write (hname, '(a, i0)') 'harm', stats%fh_harm(ihh)
    call H5Gcreate_f (g_id, trim(hname), b_id, h5_err)
    call hdf5_write_dataset_real (b_id, 'centroid', stats%fh_centroid(:,:,1:ir,ihh), err);  if (err) return
    call hdf5_write_dataset_real (b_id, 'sigma', stats%fh_sigma(:,:,1:ir,ihh), err);  if (err) return
    call hdf5_write_dataset_real (b_id, 'energy', stats%fh_energy(:,1:ir,ihh), err);  if (err) return
    call hdf5_write_dataset_real (b_id, 'power', stats%fh_power(:,1:ir,ihh), err);  if (err) return
    call hdf5_write_dataset_real (b_id, 'on_axis_intensity', stats%fh_on_axis(:,1:ir,ihh), err);  if (err) return
    call hdf5_write_dataset_real (b_id, 'emit_x', stats%fh_emit_x(:,1:ir,ihh), err);  if (err) return
    call hdf5_write_dataset_real (b_id, 'emit_y', stats%fh_emit_y(:,1:ir,ihh), err);  if (err) return
    call hdf5_write_dataset_int (b_id, 'angle_moments_valid', stats%fh_angles_valid(:,1:ir,ihh), err);  if (err) return
    call H5Gclose_f (b_id, h5_err)
  enddo
endif
call H5Gclose_f (g_id, h5_err)

! Element ends: the evaluated bunch_params_struct, unpacked into named datasets.

call H5Gcreate_f (f_id, 'element_end', g_id, h5_err)
call hdf5_write_dataset_int (g_id, 'ix_ele', stats%e_ix_ele(1:stats%iend), err);  if (err) return
call hdf5_write_dataset_real (g_id, 's', stats%e_s(1:stats%iend), err);  if (err) return

call H5Gcreate_f (g_id, 'bunch', b_id, h5_err)
call write_bp_group (b_id, stats%e_bunch(:, 1:stats%iend), err);  if (err) return
call H5Gclose_f (b_id, h5_err)

call H5Gcreate_f (g_id, 'slice', s_id, h5_err)
call write_bp_group_slices (s_id, stats%e_slice(:, :, 1:stats%iend), err);  if (err) return
call H5Gclose_f (s_id, h5_err)

! The field and bunching at each element end, so the element_end group stands alone
! when the comb keeps no per-record rows. Power/energy/intensity are totals over
! polarizations, as the per-record datasets are.

call H5Gcreate_f (g_id, 'field', s_id, h5_err)
call hdf5_write_dataset_real (s_id, 'power', stats%e_f_power(:, 1:stats%iend), err);  if (err) return
call hdf5_write_dataset_real (s_id, 'energy', stats%e_f_energy(:, 1:stats%iend), err);  if (err) return
call hdf5_write_dataset_real (s_id, 'on_axis_intensity', stats%e_f_on_axis(:, 1:stats%iend), err);  if (err) return
call hdf5_write_dataset_real (s_id, 'bunching', stats%e_bunching(:, 1:stats%iend), err);  if (err) return
call hdf5_write_dataset_real (s_id, 'bunching_phase', stats%e_bunching_phase(:, 1:stats%iend), err);  if (err) return
call H5Gclose_f (s_id, h5_err)
call H5Gclose_f (g_id, h5_err)

call H5Fclose_f (f_id, h5_err)
err_flag = .false.

!------------------------------------------------------------------------------
contains

!+
! Subroutine write_bp_group (id, rows, err)
!
! Routine to write the whole-bunch bunch_params rows as named datasets under
! group id.
!-

subroutine write_bp_group (id, rows, err)

integer(hid_t) id, mm_id
real(rp) rows(:,:)
integer im, k, h5e
logical err

!

err = .true.
call hdf5_write_dataset_real (id, 'centroid', rows(1:6, :), err);  if (err) return
call hdf5_write_dataset_real (id, 'sigma', rows(7:42, :), err);  if (err) return
call hdf5_write_dataset_real (id, 'charge_live', rows(43, :), err);  if (err) return
call hdf5_write_dataset_real (id, 'n_particle_live', rows(44, :), err);  if (err) return
call hdf5_write_dataset_real (id, 'twiss_valid', rows(45, :), err);  if (err) return
k = 45
do im = 1, 6
  call H5Gcreate_f (id, fel_stats_modes$(im), mm_id, h5e)
  call write_twiss (mm_id, rows, k, err);  if (err) return
  call H5Gclose_f (mm_id, h5e)
  k = k + 9
enddo
err = .false.

end subroutine write_bp_group

!+
! Subroutine write_twiss (id, rows, k, err)
!
! Routine to write the nine twiss datasets of one mode, reading rows(k+1:k+9, :).
!-

subroutine write_twiss (id, rows, k, err)

integer(hid_t) id
real(rp) rows(:,:)
integer k, jp
logical err

do jp = 1, 9
  call hdf5_write_dataset_real (id, trim(fel_stats_twiss$(jp)), rows(k+jp, :), err)
  if (err) return
enddo

end subroutine write_twiss

!+
! Subroutine write_bp_group_slices (id, rows, err)
!
! Routine to write the per-slice bunch_params rows as named datasets under group id.
!-

subroutine write_bp_group_slices (id, rows, err)

integer(hid_t) id, mm_id
real(rp) rows(:,:,:)
integer im, k, jp, h5e
logical err

!

err = .true.
call hdf5_write_dataset_real (id, 'centroid', rows(1:6, :, :), err);  if (err) return
call hdf5_write_dataset_real (id, 'sigma', rows(7:42, :, :), err);  if (err) return
call hdf5_write_dataset_real (id, 'charge_live', rows(43, :, :), err);  if (err) return
call hdf5_write_dataset_real (id, 'n_particle_live', rows(44, :, :), err);  if (err) return
call hdf5_write_dataset_real (id, 'twiss_valid', rows(45, :, :), err);  if (err) return
k = 45
do im = 1, 6
  call H5Gcreate_f (id, fel_stats_modes$(im), mm_id, h5e)
  do jp = 1, 9
    call hdf5_write_dataset_real (mm_id, trim(fel_stats_twiss$(jp)), rows(k+jp, :, :), err)
    if (err) return
  enddo
  call H5Gclose_f (mm_id, h5e)
  k = k + 9
enddo
err = .false.

end subroutine write_bp_group_slices

end subroutine fel_stats_write

end module fel_stats_mod
