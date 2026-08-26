!+
! Module fel_stats_mod
!
! The tracker's production statistics file <out_root>.stats.h5 (manual sec:stats).
!
! THE FILE DESCRIBES ITSELF. Every dataset goes through fel_h5_mod and carries @unit,
! @description and @axes, the last naming the coords/ datasets its dimensions run over.
! Units are fixed Bmad units (m, rad, eV, s, C, J, W) and the attributes are
! DOCUMENTATION: a reader must not scale by them.
!
! Five groups, and each answers one question.
!
!   coords/  Every axis, once. z (the one record axis), ix_ele and ele_name beside it,
!            at_element_end marking which records are element ends, and the slice axis
!            s_slice and t_slice. s_slice carries @head_direction, publishing the
!            migration invariant (higher slice index is the window head) that a
!            per-slice plot cannot be drawn without.
!   params/  Every scalar as data: the window (lambda0, window_sample, slice_spacing,
!            nbins), p0c, the charge, the species, the seed, the grid, the counts.
!            Nothing a reader needs is left in the echoed namelist.
!   beam/    slice/ holds the per-record sufficient statistics, named EXACTLY as
!            bunch_params_struct components, enough to construct one from any
!            (record, slice): centroid (nrec, nslice, 6), sigma (nrec, nslice, 6, 6) at
!            its natural rank, charge_live, n_particle_live, t, sigma_t, bunching,
!            bunching_phase, plus current and energy which every consumer would
!            otherwise re-derive. slice_twiss/ and bunch/ hold the fully evaluated Bmad
!            bunch_params, per slice and for the whole window, on the ELEMENT-END grid,
!            aligned with coords/z(at_element_end).
!   field/   total/ (always written, both polarizations of one wavelength), x/, y/ when
!            a second polarization is live, harm<h>/ per harmonic, all carrying the same
!            dataset names and their own total. No dataset's meaning depends on what
!            else the file holds.
!   meta/    Provenance, in attributes (fel_write_meta), and lattice/ beside it for
!            layout plots (fel_write_lattice).
!
! Theta moments cost FFTs, so they fill at element ends only, and angle_moments_valid
! says where (the twiss_valid pattern). Pulse-level values are POOLED downstream
! (scripts), never stored: the file stays raw.
!
! ONE RECORD AXIS. Element ends are always recorded whatever the comb (fel_comb_take),
! so an element end IS a record and the element-end arrays need no axis of their own:
! at_element_end selects them. scripts/read_stats.py is the reader everything uses, and
! scripts/bunch_params_from_stats.py reconstructs a bunch_params dict from any (record,
! slice) of the sufficient statistics. The harness checks that at element ends it
! reproduces the stored twiss.
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
use fel_h5_mod
use, intrinsic :: ieee_arithmetic

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
  integer, allocatable :: ix_ele(:)               ! (nrec) Lattice element of the record.
  integer, allocatable :: at_end(:)               ! (nrec) 1 where the record is an element
                                                  !   end. THE mask: element ends are always
                                                  !   recorded, so the element-end arrays
                                                  !   below need no axis of their own.
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
  ! Element ends: the evaluated bunch_params_struct, whole bunch and per slice. The
  ! rows align with the records at_end marks, so they carry no z or ix_ele of their
  ! own, and nothing here duplicates a per-record array.
  real(rp), allocatable :: e_bunch(:,:)           ! (n_bp, nend)
  real(rp), allocatable :: e_slice(:,:,:)         ! (n_bp, nslice, nend)
end type

!+
! Structure fel_stats_params_struct
!
! The run parameters the stats file states as DATA, so nothing downstream scrapes them
! out of the echoed namelist. Filled by the caller (fel_io_mod, which sees the whole
! run) and written into params/ verbatim.
!-

type fel_stats_params_struct
  real(rp) :: lambda0 = 0          ! Fundamental radiation wavelength [m].
  real(rp) :: slice_spacing = 0    ! Longitudinal slice spacing [m] (window_sample * lambda0).
  real(rp) :: bunch_charge = 0     ! Total live charge [C].
  real(rp) :: grid_half_width = 0  ! Transverse grid half width [m].
  real(rp) :: beta0 = 1            ! Reference beta, for the slice time axis.
  integer :: window_sample = 0     ! Slice spacing in wavelengths (Genesis's sample).
  integer :: nbins = 0             ! Beamlet size.
  integer :: grid_n_pts = 0        ! Transverse grid points per side.
  integer :: ran_seed = 0          ! The run's random seed.
  character(20) :: species = ''    ! Particle species name.
end type

! The flattened bunch_params record: 6 centroid + 36 sigma + charge_live +
! n_particle_live + twiss_valid + 6 modes x 9 params = 99.
integer, parameter :: fel_stats_n_bp$ = 99
character(*), parameter :: fel_stats_modes$(6) = [character(1):: 'x', 'y', 'z', 'a', 'b', 'c']
character(9), parameter :: fel_stats_twiss$(9) = [character(9):: &
        'beta', 'alpha', 'gamma', 'emit', 'norm_emit', 'sigma', 'sigma_p', 'eta', 'etap']
! The unit of each twiss quantity above, in the same order, for the file's @unit
! attributes. beta is in meters, the emittances in meter radians, alpha and gamma and
! the dispersions dimensionless or in meters as Bmad defines them.
character(6), parameter :: fel_stats_tunit$(9) = [character(6):: &
        'm', '1', '1/m', 'm rad', 'm rad', 'm', '1', 'm', '1']

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

allocate (stats%z(nrec), stats%ix_ele(nrec), stats%at_end(nrec))
stats%ix_ele = 0;  stats%at_end = 0
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
allocate (stats%e_bunch(fel_stats_n_bp$, nend), stats%e_slice(fel_stats_n_bp$, nslice, nend))

end subroutine fel_stats_init

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_stats_record (stats, beam, ff, z_now, ix_ele, with_angles, bdiag_arr, fpow, fonax, err_flag)
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
!   ix_ele        -- integer: Lattice element the record sits in. Written to coords/
!                      beside z, and what a layout plot joins the lattice table on.
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

subroutine fel_stats_record (stats, beam, ff, z_now, ix_ele, with_angles, bdiag_arr, fpow, fonax, err_flag)

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
integer ix_ele
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
stats%ix_ele(ir) = ix_ele
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
! Subroutine fel_stats_element_end (stats, beam, z_now, err_flag)
!
! Routine to take one element-end row: Bmad's calc_bunch_params for the whole window
! (all slices as one bunch in global window coordinates) and per slice (the Tao
! end-of-element pattern), using Bmad's own statistics machinery.
!
! The row carries the TWISS and nothing else. Everything else an element end has is
! already in the record at the same z, since an element end is always recorded
! (fel_comb_take), and this routine MARKS that record rather than copying it. The
! coincidence is a structural invariant, so it is asserted here rather than assumed.
!
! Input:
!   stats     -- fel_stats_struct: Accumulator. The current record supplies the moments.
!   beam      -- fel_beam_struct: The sliced beam, for the window chart and p0c.
!   z_now     -- real(rp): Position of the element end [m].
!
! Output:
!   stats     -- fel_stats_struct: Row stats%iend filled, at_end set on the record.
!   err_flag  -- logical: Set True on a conversion error or row overflow. False otherwise.
!-

subroutine fel_stats_element_end (stats, beam, z_now, err_flag)

type (fel_stats_struct) stats
type (fel_beam_struct), target :: beam
type (bunch_params_struct) bp
real(rp), allocatable :: cen_s(:,:), sig_s(:,:,:), w_s(:), dz_s(:), shear_s(:)
integer, allocatable :: n_s(:)
real(rp) z_now, p0_mc, p_mc
integer ie, is, nslice
logical err_flag, error
character(*), parameter :: r_name = 'fel_stats_element_end'

!

err_flag = .true.
if (stats%iend >= stats%nend) then
  call out_io (s_error$, r_name, 'MORE ELEMENT ENDS THAN THE PRECOMPUTED COUNT. INTERNAL BUG.')
  return
endif

! The invariant this row rests on: the record at this z exists, and IS this element
! end. fel_comb_take takes a row at every element end whatever the comb, so a walk
! that reordered its stats call against its element-end call would land here.

if (stats%irec < 1) then
  call out_io (s_error$, r_name, 'AN ELEMENT END WITH NO RECORD. INTERNAL BUG.')
  return
endif
if (stats%z(stats%irec) /= z_now) then
  call out_io (s_error$, r_name, 'THE LAST RECORD IS NOT AT THIS ELEMENT END. INTERNAL BUG.', &
               'RECORD \es16.8\ AGAINST ELEMENT END \es16.8\ ', &
               r_array = [stats%z(stats%irec), z_now])
  return
endif

stats%iend = stats%iend + 1
ie = stats%iend
stats%at_end(stats%irec) = 1

nslice = size(beam%slice)
allocate (cen_s(6,nslice), sig_s(6,6,nslice), w_s(nslice), dz_s(nslice), shear_s(nslice), n_s(nslice))

! ONE source of per-slice moments, the record at this z. Its two-pass weighted
! arithmetic is then the only copy of that code in the module, so nothing here can
! drift from what the per-record datasets hold.

do is = 1, nslice
  cen_s(:,is) = stats%b_centroid(:, is, stats%irec)
  sig_s(:,:,is) = reshape(stats%b_sigma(:, is, stats%irec), [6, 6])
  w_s(is) = stats%charge_live(is, stats%irec)
  n_s(is) = stats%n_particle_live(is, stats%irec)
enddo

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

! Not evaluated is NaN, not zero. A slice at the window's edge is degenerate by
! construction (no charge, collapsed modes), and a beta of zero there reads as an
! answer. twiss_valid says which case it is, and so do the numbers.

if (.not. bp%twiss_valid) then
  row(46:) = ieee_value(1.0_rp, ieee_quiet_nan)
  return
endif

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
! Subroutine fel_stats_exit_light (stats, pow, ene, bun)
!
! Routine to return the radiation the run holds at its last record: total power and
! total window energy over the live polarizations, and the mean bunching over slices.
!
! The ONE place that sums the polarizations for display, so the progress line and the
! completion block cannot come to differ about what power means. Zero before the first
! record, which only happens on a run that took none.
!
! Input:
!   stats -- fel_stats_struct: Accumulator, filled through record stats%irec.
!
! Output:
!   pow   -- real(rp): Total radiation power [W].
!   ene   -- real(rp): Total window field energy [J].
!   bun   -- real(rp): Mean |b| over slices.
!-

subroutine fel_stats_exit_light (stats, pow, ene, bun)

type (fel_stats_struct) stats
real(rp) pow, ene, bun
integer ir

!

pow = 0;  ene = 0;  bun = 0
if (stats%irec < 1) return
ir = stats%irec

pow = sum(stats%f_power(:, ir))
ene = sum(stats%f_energy(:, ir))
if (allocated(stats%f2_power)) then
  pow = pow + sum(stats%f2_power(:, ir))
  ene = ene + sum(stats%f2_energy(:, ir))
endif
bun = sum(stats%bunching(:, ir)) / stats%nslice

end subroutine fel_stats_exit_light

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_stats_write (stats, prm, file_name, err_flag)
!
! Routine to write the stats file: coords/, params/, beam/ and field/ (module header).
! meta/ and lattice/ are appended afterwards by fel_io_mod, which can see the run.
!
! Every dataset goes through fel_h5_mod, so every dataset carries @unit, @description
! and @axes. Shapes appear to h5py with the record index first: a Fortran
! (a, nslice, nrec) array reads as (nrec, nslice, a).
!
! Input:
!   stats      -- fel_stats_struct: The filled accumulator.
!   prm        -- fel_stats_params_struct: The run parameters, written to params/.
!   file_name  -- character(*): File to write.
!
! Output:
!   err_flag   -- logical: Set True on a write error. False otherwise.
!-

subroutine fel_stats_write (stats, prm, file_name, err_flag)

type (fel_stats_struct) stats
type (fel_stats_params_struct) prm
integer(hid_t) f_id, g_id, b_id, s_id
integer h5_err, ir, ie, ns, ihh, is
real(rp), allocatable :: s_slice(:), t_slice(:), cur(:,:), energy(:,:), sig_energy(:,:)
logical err_flag, err
character(*) file_name
character(12) hname
character(*), parameter :: r_name = 'fel_stats_write'

!

err_flag = .true.
err = .false.
ir = stats%irec
ie = stats%iend
ns = stats%nslice

call hdf5_open_file (file_name, 'WRITE', f_id, err);  if (err) return

call hdf5_write_attribute_string (f_id, 'file_format', 'lucifer-stats', err)
call hdf5_write_attribute_string (f_id, 'file_format_version', '2.0', err)
call hdf5_write_attribute_string (f_id, 'phase_space', 'bmad', err)
call hdf5_write_attribute_string (f_id, 'units_note', &
        'Fixed Bmad units: m, rad, eV, s, C, J, W. Every dataset carries @unit as ' // &
        'DOCUMENTATION: the values are already SI and eV, so a reader must not scale by it.', err)
if (err) return

! ------------------------------------------------------------------
! coords/: every axis, once. The record axis carries the element index and the
! element-end mask beside it, so the element-end rows below need no axis of their own.
! The element NAME is not here: it lives once in lattice/name, which ix_ele indexes.

allocate (s_slice(ns), t_slice(ns))
do is = 1, ns
  s_slice(is) = (is - 1) * prm%slice_spacing
enddo
t_slice = -s_slice / (prm%beta0 * c_light)

call H5Gcreate_f (f_id, 'coords', g_id, h5_err)
call fel_h5_real (g_id, 'z', 'm', 'Path length along the line at each record.', '', stats%z(1:ir), err)
call fel_h5_int (g_id, 'ix_ele', '1', &
      'Lattice element holding each record. Indexes the lattice/ table.', 'z', stats%ix_ele(1:ir), err)
call fel_h5_flag (g_id, 'at_element_end', &
      'One where the record is an element end. Selects the beam/bunch and ' // &
      'beam/slice_twiss rows, in order.', 'z', stats%at_end(1:ir), err)
call fel_h5_real (g_id, 's_slice', 'm', &
      'Slice position in the time window, slice 1 at zero.', '', s_slice, err)
call fel_h5_real (g_id, 't_slice', 's', &
      'Arrival-time offset of each slice, negative toward the window head.', '', t_slice, err)
if (err) return

! The head convention, which no per-slice plot can be drawn without: the migration
! invariant z_global = z_local + beta*(islice-1)*spacing puts the HEAD at the high
! index, so s_slice increases toward the head.

call H5LTset_attribute_string_f (g_id, 's_slice', 'head_direction', '+index', h5_err)
err = err .or. (h5_err < 0)
call H5Gclose_f (g_id, h5_err)
if (err) return

! ------------------------------------------------------------------
! params/: every scalar as data.

call H5Gcreate_f (f_id, 'params', g_id, h5_err)
call fel_h5_real (g_id, 'lambda0', 'm', 'Fundamental radiation wavelength.', '', prm%lambda0, err)
call fel_h5_int (g_id, 'window_sample', '1', &
      'Slice spacing in wavelengths (Genesis''s sample).', '', prm%window_sample, err)
call fel_h5_real (g_id, 'slice_spacing', 'm', &
      'Longitudinal slice spacing, window_sample * lambda0.', '', prm%slice_spacing, err)
call fel_h5_int (g_id, 'nbins', '1', 'Beamlet size of the quiet loading.', '', prm%nbins, err)
call fel_h5_real (g_id, 'p0c', 'eV', 'Reference momentum times c.', '', stats%p0c, err)
call fel_h5_real (g_id, 'bunch_charge', 'C', 'Total live charge of the window.', '', prm%bunch_charge, err)
call fel_h5_int (g_id, 'ran_seed', '1', 'Random seed of the run.', '', prm%ran_seed, err)
call fel_h5_int (g_id, 'grid_n_pts', '1', 'Transverse grid points per side.', '', prm%grid_n_pts, err)
call fel_h5_real (g_id, 'grid_half_width', 'm', 'Transverse grid half width.', '', prm%grid_half_width, err)
call fel_h5_int (g_id, 'n_record', '1', 'Records taken.', '', ir, err)
call fel_h5_int (g_id, 'n_element_end', '1', &
      'Element ends taken. Equals the count of coords/at_element_end.', '', ie, err)
call fel_h5_int (g_id, 'n_slice', '1', 'Slices in the time window.', '', ns, err)
call fel_h5_str (g_id, 'species', 'Particle species.', '', [prm%species], err)
call H5Gclose_f (g_id, h5_err)
if (err) return

! ------------------------------------------------------------------
! beam/slice/: the per-record sufficient statistics, bunch_params_struct names, plus
! the two derived quantities every consumer would otherwise re-derive.

allocate (cur(ns, ir), energy(ns, ir), sig_energy(ns, ir))
cur = c_light * stats%charge_live(:,1:ir) / prm%slice_spacing
energy = sqrt((stats%p0c * (1 + stats%b_centroid(6,:,1:ir)))**2 + m_electron**2)
! sigma_E = beta * p0c * sigma_pz, since dE/dp = p/E = beta. The factor is 1 to a few
! parts in 1e9 at an X-ray FEL's energy and it is not 1 for everything this tracker
! could be pointed at, so it is here rather than dropped.
sig_energy = stats%p0c * sqrt(max(0.0_rp, stats%b_sigma(36,:,1:ir))) * &
             stats%p0c * (1 + stats%b_centroid(6,:,1:ir)) / energy

call H5Gcreate_f (f_id, 'beam', b_id, h5_err)
call H5Gcreate_f (b_id, 'slice', g_id, h5_err)
call fel_h5_real (g_id, 'centroid', 'm,1,m,1,m,1', &
      'Weighted centroid, Bmad phase space (x, px/p0, y, py/p0, z, pz).', &
      'z,s_slice', stats%b_centroid(:,:,1:ir), err)
call fel_h5_real (g_id, 'sigma', 'm,1,m,1,m,1 squared', &
      'Second moments about the centroid, the 6x6 at its natural rank.', &
      'z,s_slice', reshape(stats%b_sigma(:,:,1:ir), [6, 6, ns, ir]), err)
call fel_h5_real (g_id, 'charge_live', 'C', 'Live charge of the slice.', &
      'z,s_slice', stats%charge_live(:,1:ir), err)
call fel_h5_int (g_id, 'n_particle_live', '1', 'Live macroparticles in the slice.', &
      'z,s_slice', stats%n_particle_live(:,1:ir), err)
call fel_h5_real (g_id, 't', 's', &
      'Mean arrival time of the slice at this plane, relative to the reference.', &
      'z,s_slice', stats%t(:,1:ir), err)
call fel_h5_real (g_id, 'sigma_t', 's', 'Rms arrival-time spread within the slice.', &
      'z,s_slice', stats%sigma_t(:,1:ir), err)
call fel_h5_real (g_id, 'bunching', '1', &
      'Bunching |b| at the fundamental, charge weighted.', 'z,s_slice', stats%bunching(:,1:ir), err)
call fel_h5_real (g_id, 'bunching_phase', 'rad', &
      'Bunching phase arg(b), carrying the run''s reference phase.', &
      'z,s_slice', stats%bunching_phase(:,1:ir), err)
call fel_h5_real (g_id, 'current', 'A', &
      'Slice current, c * charge_live / slice_spacing.', 'z,s_slice', cur, err)
call fel_h5_real (g_id, 'energy', 'eV', &
      'Mean total energy of the slice, Bmad''s convention (energy is eV, never gamma).', &
      'z,s_slice', energy, err)
call fel_h5_real (g_id, 'sigma_energy', 'eV', 'Rms energy spread of the slice.', &
      'z,s_slice', sig_energy, err)
call H5Gclose_f (g_id, h5_err)
if (err) return

! beam/slice_twiss/ and beam/bunch/: the evaluated bunch_params, on the ELEMENT-END
! grid. They align with the records coords/at_element_end marks, in order.

call H5Gcreate_f (b_id, 'slice_twiss', g_id, h5_err)
call write_bp_slices (g_id, stats%e_slice(:, :, 1:ie), err)
call H5Gclose_f (g_id, h5_err)
if (err) return

call H5Gcreate_f (b_id, 'bunch', g_id, h5_err)
call write_bp_bunch (g_id, stats%e_bunch(:, 1:ie), err)
call H5Gclose_f (g_id, h5_err)
call H5Gclose_f (b_id, h5_err)
if (err) return

! ------------------------------------------------------------------
! field/: total/ always, then one group per component and per harmonic, all carrying
! the same dataset names. Nothing here changes meaning when a sibling appears.

call H5Gcreate_f (f_id, 'field', b_id, h5_err)

if (allocated(stats%f2_power)) then
  call write_field_total (b_id, stats%f_power(:,1:ir) + stats%f2_power(:,1:ir), &
                          stats%f_energy(:,1:ir) + stats%f2_energy(:,1:ir), &
                          stats%f_on_axis(:,1:ir) + stats%f2_on_axis(:,1:ir), 1, err)
else
  call write_field_total (b_id, stats%f_power(:,1:ir), stats%f_energy(:,1:ir), &
                          stats%f_on_axis(:,1:ir), 1, err)
endif
if (err) return

call H5Gcreate_f (b_id, 'x', g_id, h5_err)
call write_field_component (g_id, stats%f_centroid(:,:,1:ir), stats%f_sigma(:,:,1:ir), &
        stats%f_power(:,1:ir), stats%f_energy(:,1:ir), stats%f_on_axis(:,1:ir), &
        stats%f_emit_x(:,1:ir), stats%f_emit_y(:,1:ir), stats%f_angles_valid(:,1:ir), err)
call H5Gclose_f (g_id, h5_err)
if (err) return

if (allocated(stats%f2_power)) then
  call H5Gcreate_f (b_id, 'y', g_id, h5_err)
  call write_field_component (g_id, stats%f2_centroid(:,:,1:ir), stats%f2_sigma(:,:,1:ir), &
          stats%f2_power(:,1:ir), stats%f2_energy(:,1:ir), stats%f2_on_axis(:,1:ir), &
          stats%f2_emit_x(:,1:ir), stats%f2_emit_y(:,1:ir), stats%f_angles_valid(:,1:ir), err)
  call H5Gclose_f (g_id, h5_err)
  if (err) return
endif

! Harmonics: a group each, with its own total, since a detector separates colors and
! nothing may sum across wavelengths.

if (allocated(stats%fh_power)) then
  do ihh = 1, size(stats%fh_harm)
    write (hname, '(a, i0)') 'harm', stats%fh_harm(ihh)
    call H5Gcreate_f (b_id, trim(hname), s_id, h5_err)
    call H5LTset_attribute_int_f (b_id, trim(hname), 'harmonic', [stats%fh_harm(ihh)], 1_size_t, h5_err)
    call write_field_total (s_id, stats%fh_power(:,1:ir,ihh), stats%fh_energy(:,1:ir,ihh), &
                            stats%fh_on_axis(:,1:ir,ihh), stats%fh_harm(ihh), err)
    call H5Gcreate_f (s_id, 'x', g_id, h5_err)
    call write_field_component (g_id, stats%fh_centroid(:,:,1:ir,ihh), stats%fh_sigma(:,:,1:ir,ihh), &
            stats%fh_power(:,1:ir,ihh), stats%fh_energy(:,1:ir,ihh), stats%fh_on_axis(:,1:ir,ihh), &
            stats%fh_emit_x(:,1:ir,ihh), stats%fh_emit_y(:,1:ir,ihh), stats%fh_angles_valid(:,1:ir,ihh), err)
    call H5Gclose_f (g_id, h5_err)
    call H5Gclose_f (s_id, h5_err)
    if (err) return
  enddo
endif

call H5Gclose_f (b_id, h5_err)

call H5Fclose_f (f_id, h5_err)
err_flag = .false.

!------------------------------------------------------------------------------
contains

!+
! Subroutine write_field_total (id, pow, ene, onax, harm, err)
!
! Routine to write the total/ group of one wavelength: the sum over the live
! polarizations. ALWAYS written, whether one polarization is live or two, so that no
! reader has to ask what else the file holds before it knows what power means.
!-

subroutine write_field_total (id, pow, ene, onax, harm, err)

integer(hid_t) id, tt_id
integer harm, h5e
real(rp) pow(:,:), ene(:,:), onax(:,:)
logical err
character(60) note

!

write (note, '(a, i0, a)') ' Summed over the live polarizations of harmonic ', harm, '.'

call H5Gcreate_f (id, 'total', tt_id, h5e)
call fel_h5_real (tt_id, 'power', 'W', 'Radiation power of the slice.' // trim(note), &
      'z,s_slice', pow, err)
call fel_h5_real (tt_id, 'energy', 'J', 'Field energy of the slice.' // trim(note), &
      'z,s_slice', ene, err)
call fel_h5_real (tt_id, 'on_axis_intensity', 'W/m^2', &
      'Intensity at the grid center.' // trim(note), 'z,s_slice', onax, err)
call H5Gclose_f (tt_id, h5e)

end subroutine write_field_total

!+
! Subroutine write_field_component (id, cen, sig, pow, ene, onax, ex, ey, valid, err)
!
! Routine to write one polarization component: the full wavefront_params set under the
! same dataset names every component and harmonic uses.
!-

subroutine write_field_component (id, cen, sig, pow, ene, onax, ex, ey, valid, err)

integer(hid_t) id
real(rp) cen(:,:,:), sig(:,:,:), pow(:,:), ene(:,:), onax(:,:), ex(:,:), ey(:,:)
integer valid(:,:)
logical err

!

call fel_h5_real (id, 'centroid', 'm,rad,m,rad', &
      'Intensity-weighted (x, theta_x, y, theta_y). The theta entries are NaN where ' // &
      'the angle moments were not computed.', 'z,s_slice', cen, err)
call fel_h5_real (id, 'sigma', 'm,rad,m,rad squared', &
      'Wigner second moments, the 4x4 at its natural rank. The theta rows are NaN ' // &
      'where the angle moments were not computed.', &
      'z,s_slice', reshape(sig, [4, 4, size(sig,2), size(sig,3)]), err)
call fel_h5_real (id, 'power', 'W', 'Radiation power of this component.', 'z,s_slice', pow, err)
call fel_h5_real (id, 'energy', 'J', 'Field energy of this component.', 'z,s_slice', ene, err)
call fel_h5_real (id, 'on_axis_intensity', 'W/m^2', &
      'Intensity of this component at the grid center.', 'z,s_slice', onax, err)
call fel_h5_real (id, 'emit_x', 'm rad', &
      'sqrt(det sigma_x-plane) = M^2 lambda / 4 pi. NaN without angle moments.', &
      'z,s_slice', ex, err)
call fel_h5_real (id, 'emit_y', 'm rad', &
      'sqrt(det sigma_y-plane) = M^2 lambda / 4 pi. NaN without angle moments.', &
      'z,s_slice', ey, err)
call fel_h5_flag (id, 'angle_moments_valid', &
      'One where the theta moments were computed (they cost three FFTs, so element ' // &
      'ends only).', 'z,s_slice', valid, err)

end subroutine write_field_component

!+
! Subroutine write_bp_bunch (id, rows, err)
!
! Routine to write the whole-window bunch_params rows: the packed 99-number row
! unpacked into named datasets, twiss modes in their own groups.
!-

subroutine write_bp_bunch (id, rows, err)

integer(hid_t) id, mm_id
real(rp) rows(:,:)
integer im, k, jp, h5e, nr
logical err

!

nr = size(rows, 2)
call fel_h5_real (id, 'centroid', 'm,1,m,1,m,1', &
      'Weighted centroid of the whole window, Bmad phase space.', &
      'element_end', rows(1:6, :), err)
call fel_h5_real (id, 'sigma', 'm,1,m,1,m,1 squared', &
      'Second moments of the whole window, the 6x6 at its natural rank.', &
      'element_end', reshape(rows(7:42, :), [6, 6, nr]), err)
call fel_h5_real (id, 'charge_live', 'C', 'Live charge of the whole window.', &
      'element_end', rows(43, :), err)
call fel_h5_int (id, 'n_particle_live', '1', 'Live macroparticles in the whole window.', &
      'element_end', nint(rows(44, :)), err)
call fel_h5_flag (id, 'twiss_valid', &
      'One where Bmad evaluated the twiss (it needs six live particles).', &
      'element_end', nint(rows(45, :)), err)
k = 45
do im = 1, 6
  call H5Gcreate_f (id, fel_stats_modes$(im), mm_id, h5e)
  do jp = 1, 9
    call fel_h5_real (mm_id, trim(fel_stats_twiss$(jp)), trim(fel_stats_tunit$(jp)), &
          trim(fel_stats_modes$(im)) // '-mode ' // trim(fel_stats_twiss$(jp)) // &
          '. NaN where twiss_valid is zero.', 'element_end', rows(k+jp, :), err)
  enddo
  call H5Gclose_f (mm_id, h5e)
  k = k + 9
enddo

end subroutine write_bp_bunch

!+
! Subroutine write_bp_slices (id, rows, err)
!
! Routine to write the per-slice bunch_params rows, the same names as write_bp_bunch
! with the slice axis added.
!-

subroutine write_bp_slices (id, rows, err)

integer(hid_t) id, mm_id
real(rp) rows(:,:,:)
integer im, k, jp, h5e, nsl, nr
logical err

!

nsl = size(rows, 2);  nr = size(rows, 3)
call fel_h5_real (id, 'centroid', 'm,1,m,1,m,1', &
      'Weighted centroid per slice, Bmad phase space.', &
      'element_end,s_slice', rows(1:6, :, :), err)
call fel_h5_real (id, 'sigma', 'm,1,m,1,m,1 squared', &
      'Second moments per slice, the 6x6 at its natural rank.', &
      'element_end,s_slice', reshape(rows(7:42, :, :), [6, 6, nsl, nr]), err)
call fel_h5_real (id, 'charge_live', 'C', 'Live charge of the slice.', &
      'element_end,s_slice', rows(43, :, :), err)
call fel_h5_int (id, 'n_particle_live', '1', 'Live macroparticles in the slice.', &
      'element_end,s_slice', nint(rows(44, :, :)), err)
call fel_h5_flag (id, 'twiss_valid', &
      'One where Bmad evaluated the twiss (it needs six live particles).', &
      'element_end,s_slice', nint(rows(45, :, :)), err)
k = 45
do im = 1, 6
  call H5Gcreate_f (id, fel_stats_modes$(im), mm_id, h5e)
  do jp = 1, 9
    call fel_h5_real (mm_id, trim(fel_stats_twiss$(jp)), trim(fel_stats_tunit$(jp)), &
          trim(fel_stats_modes$(im)) // '-mode ' // trim(fel_stats_twiss$(jp)) // &
          ' per slice. NaN where twiss_valid is zero.', 'element_end,s_slice', rows(k+jp, :, :), err)
  enddo
  call H5Gclose_f (mm_id, h5e)
  k = k + 9
enddo

end subroutine write_bp_slices

end subroutine fel_stats_write

end module fel_stats_mod
