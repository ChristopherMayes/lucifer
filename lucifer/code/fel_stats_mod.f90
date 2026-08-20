!+
! Module fel_stats_mod
!
! The tracker's production statistics file <out_root>.stats.h5 (manual sec:stats).
!
! Fixed Bmad units everywhere -- m, rad, eV, s, C, J, W -- with units attributes written
! as documentation only, never load-bearing. Per-record datasets are (nrec, nslice)
! arrays as read by h5py (Genesis4-style per-slice export for visualization), with beam
! datasets named EXACTLY as bunch_params_struct components, sufficient to construct one
! from any (record, slice): centroid (nrec, nslice, 6), sigma (nrec, nslice, 36 -- the
! 6x6 flattened; symmetric, so the flattening order cannot mislead), charge_live,
! charge_tot, n_particle_live, n_particle_tot, t, sigma_t, plus z(nrec) and s == z.
! The field side stores wavefront_params_struct components per slice: centroid
! (nrec, nslice, 4), sigma (nrec, nslice, 16), energy, power, on_axis_intensity,
! emit_x, emit_y, angle_moments_valid (0/1; theta moments cost FFTs, so they fill at
! element ends only -- the twiss_valid pattern). Pulse-level values are POOLED
! downstream (scripts), never stored: the file stays raw.
!
! At ELEMENT ENDS the fully evaluated Bmad bunch_params_struct lands under
! element_end/: per whole-window bunch AND per slice, via Bmad's own calc_bunch_params
! (the Tao end-of-element pattern) -- twiss groups x/y/z/a/b/c each carrying beta,
! alpha, gamma, emit, norm_emit, sigma, sigma_p, eta, etap; plus centroid, sigma,
! charge_live, n_particle_live, ix_ele, s, twiss_valid. scripts/
! bunch_params_from_stats.py reconstructs a bunch_params dict from any (record, slice)
! of the per-record sufficient statistics; the harness checks that at element ends it
! reproduces these stored values.
!
! Beam moments are computed two-pass (mean first, then centered second moments -- the
! FINDINGS 4.8 variance lesson), weighted by macroparticle charge, parallel over
! slices with fixed-order results (each slice's sums are its own).
!
! t and sigma_t derive from the stored chart: t = z_now/c - <z>/(beta0 c) and
! sigma_t = sigma_z/(beta0 c) with beta0 the reference beta -- exact for the reference
! particle, and documented rather than hidden inside a convention.
!-

module fel_stats_mod

use bmad
use beam_utils, only: calc_bunch_params, calc_emittances_and_twiss_from_sigma_matrix
use fel_beam_mod
use wavefront_mod
use fel_track_mod, only: fel_slip_struct, fel_field_index, fel_field_diag
use hdf5_interface

implicit none

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
  ! the whole story; with two, they carry the X component, the f2_* arrays carry Y
  ! (written as a field/y/ group), and the power/energy/on_axis datasets are written
  ! as TOTALS -- single-polarization files are unchanged.
  real(rp), allocatable :: f_centroid(:,:,:)      ! (4, nslice, nrec)
  real(rp), allocatable :: f_sigma(:,:,:)         ! (16, nslice, nrec)
  real(rp), allocatable :: f_energy(:,:), f_power(:,:), f_on_axis(:,:)
  real(rp), allocatable :: f_emit_x(:,:), f_emit_y(:,:)
  integer, allocatable :: f_angles_valid(:,:)     ! 0/1
  real(rp), allocatable :: f2_centroid(:,:,:), f2_sigma(:,:,:)
  real(rp), allocatable :: f2_energy(:,:), f2_power(:,:), f2_on_axis(:,:)
  real(rp), allocatable :: f2_emit_x(:,:), f2_emit_y(:,:)
  ! Element ends: the evaluated bunch_params_struct, whole bunch and per slice.
  integer, allocatable :: e_ix_ele(:)             ! (nend)
  real(rp), allocatable :: e_s(:)                 ! (nend) [m]
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
! Subroutine fel_stats_init (stats, nrec, nend, nslice)
!
! Routine to size the accumulation arrays. nrec and nend are exact counts precomputed
! from the lattice walk (records = 1 + sum of undulator steps + one per interlude).
!-

subroutine fel_stats_init (stats, nrec, nend, nslice, p0c, two_pol)

type (fel_stats_struct) stats
real(rp) p0c
integer nrec, nend, nslice
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
allocate (stats%e_ix_ele(nend), stats%e_s(nend))
allocate (stats%e_bunch(fel_stats_n_bp$, nend), stats%e_slice(fel_stats_n_bp$, nslice, nend))

end subroutine fel_stats_init

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_stats_record (stats, beam, wf, slip, z_now, with_angles, err_flag)
!
! Routine to take one per-record row: beam sufficient statistics and wavefront params
! for every slice. with_angles fills the field theta moments (element ends).
!-

subroutine fel_stats_record (stats, beam, wf, slip, z_now, with_angles, bdiag_arr, fpow, fonax, err_flag)

type (fel_stats_struct) stats
type (fel_beam_struct), target :: beam
type (wavefront_struct), target :: wf
type (fel_slip_struct) slip
type (fel_slice_struct), pointer :: sl
type (wavefront_params_struct) pms
type (fel_slice_diag_struct) bdiag_arr(:)   ! OUT: the diag instrument, one per slice.
real(rp) fpow(:), fonax(:)                  ! OUT: fel_field_diag's power/on-axis per slice.
real(rp) z_now
logical with_angles, err_flag

real(rp) w, wsum, mean(6), cen(6), sig(6,6), v(6), beta0, ks
integer ir, is, ip, i, j, nslice
logical err, any_err
character(*), parameter :: r_name = 'fel_stats_record'

!

err_flag = .true.
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
! every slice -- the diag writer then only prints. Each slice's arithmetic is the
! identical serial code, so diag.txt is bit-for-bit what it always was; what changed
! is that the formerly SERIAL per-record diag sweeps now ride this parallel loop.

!$OMP parallel do private(sl, w, wsum, mean, cen, sig, v, ip, i, j, pms, err) &
!$OMP&   reduction(.or.: any_err)
do is = 1, nslice
  sl => beam%slice(is)
  call fel_field_diag (wf, fel_field_index(slip, is, nslice), fpow(is), fonax(is))
  call fel_slice_diag (beam, sl, ks, bdiag_arr(is))

  ! Two-pass weighted moments (the FINDINGS 4.8 variance lesson).

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

  ! Bunching from the diag instrument just evaluated (|b| at the fundamental; the
  ! phase carries the common phi0, Genesis's own convention).

  stats%bunching(is, ir) = bdiag_arr(is)%bunching
  stats%bunching_phase(is, ir) = bdiag_arr(is)%bunching_phase

  ! The field slice this beam slice couples to, unrotated exactly as the dumps are.

  call wavefront_params_of_plane (wf%Ex(:,:,fel_field_index(slip, is, nslice)), wf%dx, &
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
    call wavefront_params_of_plane (wf%Ey(:,:,fel_field_index(slip, is, nslice)), wf%dx, &
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
! (all slices as one bunch in global window coordinates) and per slice -- the Tao
! end-of-element pattern, using Bmad's own statistics machinery.
!-

subroutine fel_stats_element_end (stats, beam, ele, z_now, err_flag)

type (fel_stats_struct) stats
type (fel_beam_struct), target :: beam
type (ele_struct) ele
type (bunch_struct) bunch
type (bunch_params_struct) bp
real(rp), allocatable :: beta0(:)
real(rp) z_now
integer ie, is
logical err_flag, err, error
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

call fel_concat_slices (beam, ele, bunch, beta0, err);  if (err) return
call calc_bunch_params (bunch, bp, error)
call pack_bp (bp, stats%e_bunch(:, ie))

! Per slice: the twiss evaluation is Bmad's own (calc_emittances_and_twiss_from_
! sigma_matrix, the identical code path calc_bunch_params ends in), but fed from the
! moments the CURRENT record already computed -- an element end always coincides with
! its last record -- instead of re-summing every particle through a bunch conversion.
! Measured: this is what moved the element-end cost from 7% of the demo run into the
! noise (2208 conversions + re-summations retired per run).

do is = 1, size(beam%slice)
  bp = bunch_params_struct()
  bp%centroid%vec = stats%b_centroid(:, is, stats%irec)
  bp%centroid%p0c = stats%p0c
  bp%centroid%species = electron$
  bp%sigma = reshape(stats%b_sigma(:, is, stats%irec), [6, 6])
  bp%charge_live = stats%charge_live(is, stats%irec)
  bp%n_particle_live = stats%n_particle_live(is, stats%irec)
  bp%n_particle_tot = bp%n_particle_live
  if (bp%n_particle_live >= 6) then
    call calc_emittances_and_twiss_from_sigma_matrix (bp%sigma, bp, error, .false.)
  endif
  call pack_bp (bp, stats%e_slice(:, is, ie))
enddo

err_flag = .false.

!------------------------------------------------------------------------------
contains

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
! Subroutine fel_stats_write (stats, file_name, err_flag)
!
! Routine to write the stats file. Datasets appear to h5py with the record index first:
! Fortran (a, nslice, nrec) reads as (nrec, nslice, a). Units attributes are written as
! documentation; the units are FIXED and the attributes are never load-bearing.
!-

subroutine fel_stats_write (stats, file_name, err_flag)

type (fel_stats_struct) stats
integer(hid_t) f_id, g_id, b_id, s_id
integer h5_err, ir
logical err_flag, err
character(*) file_name
character(*), parameter :: r_name = 'fel_stats_write'

!

err_flag = .true.
ir = stats%irec

call hdf5_open_file (file_name, 'WRITE', f_id, err);  if (err) return

call hdf5_write_dataset_real (f_id, 'z', stats%z(1:ir), err);  if (err) return
call hdf5_write_dataset_real (f_id, 'p0c', [stats%p0c], err);  if (err) return

! One root attribute documents the fixed units; it is never load-bearing.

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
  ! Two live polarizations: power/energy/intensity are TOTALS; the x-component params
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
call H5Gclose_f (g_id, h5_err)

call H5Fclose_f (f_id, h5_err)
err_flag = .false.

!------------------------------------------------------------------------------
contains

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
