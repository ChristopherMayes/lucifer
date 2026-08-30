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
!            the slice axis with ct_slice, t_slice and z_slice on it, all three carrying
!            @head_direction, publishing the
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
  ! Per-coordinate extremes RELATIVE TO THE CENTROID, bunch_params_struct's rel_max and
  ! rel_min: six phase-space entries plus time. Order statistics, which no moment can
  ! reconstruct, kept for envelope plots. NaN for an empty slice.
  real(rp), allocatable :: b_rel_max(:,:,:)       ! (7, nslice, nrec)
  real(rp), allocatable :: b_rel_min(:,:,:)       ! (7, nslice, nrec)
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
! The few run scalars the stats writer itself needs, filled by the caller (fel_io_mod,
! which sees the whole run) and written into run/. The INPUT tree params/ is written by
! fel_io_mod directly, one subgroup per honored struct, since only it sees the run.
!-

type fel_stats_params_struct
  real(rp) :: slice_spacing = 0    ! Longitudinal slice spacing [m] (window_sample * lambda0).
  real(rp) :: beta0 = 1            ! Reference beta, for the slice z variable.
  character(20) :: species = ''    ! Particle species name.
end type

! The flattened bunch_params record: 6 centroid + 36 sigma + charge_live +
! n_particle_live + twiss_valid + 6 modes x 9 params = 99.
integer, parameter :: fel_stats_n_bp$ = 99
! The PACKED order of a bunch_params row's twiss blocks. The first three are the
! projected planes and the last three the normal modes, which is why the file puts them
! on two axes: they are two decompositions of one beam, and an eigen-emittance is not a
! projected emittance.

character(*), parameter :: fel_stats_modes$(6) = [character(1):: 'x', 'y', 'z', 'a', 'b', 'c']
character(*), parameter :: fel_stats_planes$(3) = [character(1):: 'x', 'y', 'z']
character(*), parameter :: fel_stats_norm$(3) = [character(1):: 'a', 'b', 'c']

! The label axes of a vector and of a matrix, written into coords/ so that no reader has
! to guess a trailing dimension from its length, with the unit list that goes with each.

character(5), parameter :: fel_stats_bmad$(6) = [character(5):: 'x', 'px', 'y', 'py', 'z', 'pz']
character(7), parameter :: fel_stats_wf$(4) = [character(7):: 'x', 'theta_x', 'y', 'theta_y']
! The polarization component names, in the order the file writes them.
character(*), parameter :: fel_stats_pol$(2) = [character(1):: 'x', 'y']
character(*), parameter :: fel_stats_bmad_unit$ = 'm,1,m,1,m,1'
character(*), parameter :: fel_stats_wf_unit$ = 'm,rad,m,rad'
! The same two unit lists PER ENTRY, written into coords/ beside the labels they belong
! to. A unit is a property of the axis, so a dataset over that axis says only which axis
! carries its units and to what power.
character(3), parameter :: fel_stats_bmad_u$(6) = [character(3):: 'm', '1', 'm', '1', 'm', '1']
character(3), parameter :: fel_stats_wf_u$(4) = [character(3):: 'm', 'rad', 'm', 'rad']
! Every value @kind takes, enumerated on the root so a reader that sorts groups by kind
! does not have to discover the vocabulary by inspection. An ARRAY: shape expresses
! arity, and this holds a list.
character(10), parameter :: fel_stats_kinds$(14) = [character(10):: 'axis', 'input', &
        'run', 'beam', 'per_slice', 'projected', 'twiss', 'modes', 'field', &
        'component', 'derived', 'harmonic', 'table', 'provenance']
! The label axis of the per-slice extremes: the six phase-space names plus time, since
! bunch_params_struct's rel_max and rel_min are seven-vectors.
character(2), parameter :: fel_stats_bmad_t$(7) = [character(2):: 'x', 'px', 'y', 'py', 'z', 'pz', 't']
character(3), parameter :: fel_stats_bmad_t_u$(7) = [character(3):: 'm', '1', 'm', '1', 'm', '1', 's']
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
allocate (stats%b_rel_max(7, nslice, nrec), stats%b_rel_min(7, nslice, nrec))
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

real(rp) w, wsum, mean(6), cen(6), sig(6,6), v(6), vmin(6), vmax(6), beta0, ks
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

!$OMP parallel do private(sl, w, wsum, mean, cen, sig, v, vmin, vmax, ip, i, j, io, pms, err) &
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
  vmin = huge(1.0_rp);  vmax = -huge(1.0_rp)
  do ip = 1, sl%n
    v = [sl%x(ip), sl%px(ip), sl%y(ip), sl%py(ip), sl%z(ip), sl%pz(ip)] - mean
    vmin = min(vmin, v);  vmax = max(vmax, v)
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

  ! The per-coordinate extremes, relative to the centroid just computed: order
  ! statistics that ride the sweep for free. The time entry maps through
  ! t - <t> = -(z - <z>)/(beta0 c), so the LARGEST time offset is the SMALLEST z one.

  if (sl%n > 0) then
    stats%b_rel_max(1:6, is, ir) = vmax
    stats%b_rel_min(1:6, is, ir) = vmin
    stats%b_rel_max(7, is, ir) = -vmin(5) / (beta0 * c_light)
    stats%b_rel_min(7, is, ir) = -vmax(5) / (beta0 * c_light)
  else
    stats%b_rel_max(:, is, ir) = ieee_value(1.0_rp, ieee_quiet_nan)
    stats%b_rel_min(:, is, ir) = ieee_value(1.0_rp, ieee_quiet_nan)
  endif

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
! Routine to write the stats file: the identity, coords/, run/, beam/ and field/.
! params/ (the input tree), meta/, lattice/ and coords/ele are appended afterwards by
! fel_io_mod, which can see the run.
!
! Every dataset goes through fel_h5_mod, so every dataset carries @unit, @long_name,
! @description and @axes, and EVERY NAME IN @axes RESOLVES TO A coords/ DATASET,
! trailing label axes included. Shapes appear to h5py with the record index first: a
! Fortran (a, nslice, nrec) array reads as (nrec, nslice, a).
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
integer h5_err, ir, ie, ns, ihh, is, ip, ne
integer, allocatable :: rec(:), e_ix(:)
integer, allocatable :: sl(:)
real(rp), allocatable :: ct_slice(:), t_slice(:), z_slice(:), e_s(:)
real(rp), allocatable :: cur(:,:), energy(:,:), sig_energy(:,:)
logical, allocatable :: mask(:)
logical err_flag, err
character(*) file_name
logical merr
character(12) hname
character(24) writer_str
character(*), parameter :: r_name = 'fel_stats_write'

!

err_flag = .true.
err = .false.
ir = stats%irec
ie = stats%iend
ns = stats%nslice

call hdf5_open_file (file_name, 'WRITE', f_id, err);  if (err) return

! The identity (bmad-stats R1 to R5): the FORMAT is bmad-stats with the fel extension,
! the planned reset from lucifer-stats 2.x, and the writer is this program, whose only
! version is the Bmad it was built against.

call hdf5_write_attribute_string (f_id, 'file_format', 'bmad-stats', err)
call hdf5_write_attribute_string (f_id, 'file_format_version', '1.0', err)
write (writer_str, '(a, i0)') 'lucifer bmad-', bmad_inc_version$
call hdf5_write_attribute_string (f_id, 'writer', trim(writer_str), err)
call hdf5_write_attribute_string_rank1 (f_id, 'extensions', [character(3):: 'fel'], merr)
err = err .or. merr
call hdf5_write_attribute_string (f_id, 'fel_version', '1.0', err)
call hdf5_write_attribute_string (f_id, 'phase_space', 'bmad', err)
call hdf5_write_attribute_string_rank1 (f_id, 'kinds', fel_stats_kinds$, merr)
err = err .or. merr
call hdf5_write_attribute_string (f_id, 'units_note', &
        'Fixed Bmad units: m, rad, eV, s, C, J, W. Every dataset carries @unit as ' // &
        'DOCUMENTATION: the values are already SI and eV, so a reader must not scale ' // &
        'by it. Every name in @axes is a coords/ dataset. A one-byte integer dataset ' // &
        'with @unit = 1 and @dtype_hint = bool is this format''s boolean, HDF5 ' // &
        'having none. Every group carries @kind, from the root''s @kinds list. An ' // &
        'attribute holding ONE value is a scalar and one holding a list is an array, ' // &
        'length one included. A coordinate VARIABLE may repeat: coords/s does, wherever ' // &
        'a zero-length element applies a wake kick, which is why the record number and ' // &
        'not s is the axis.', err)
if (err) return

! ------------------------------------------------------------------
! coords/: every axis, once, each saying which axis it IS. The record NUMBER is the
! axis, not z: z repeats wherever two records land on one plane, and an index that
! repeats answers a selection silently wrong. z rides along as a variable on it.

allocate (mask(ir), rec(ir), sl(ns), ct_slice(ns), t_slice(ns), z_slice(ns))
do ip = 1, ir
  rec(ip) = ip - 1
enddo
mask = (stats%at_end(1:ir) == 1)
ne = count(mask)
allocate (e_ix(ne), e_s(ne))
e_ix = pack(stats%ix_ele(1:ir), mask)
e_s = pack(stats%z(1:ir), mask)

! The slice grid is uniform in TIME, and that is not a convention. fel_concat_slices
! holds z_global = z_local + beta_j*(islice-1)*spacing with the PARTICLE's beta, and
! z = -beta*c*(t - t_ref), so at a grid point t - t_ref = -(islice-1)*spacing/c and the
! beta CANCELS. fel_apply_slippage accumulates in radiation wavelengths and rotates the
! record one slice at window_sample, and the slippage rate (1+aw^2)/(2 gamma^2 lambda)
! is exactly one wavelength per undulator period, so one slice is exactly
! window_sample wavelengths of slippage and the rotation needs no interpolation. That
! is what makes the light-travel distance the field's own coordinate.
!
! So ct_slice and t_slice are EXACT and carry no beta. z_slice is Bmad's z at the
! REFERENCE beta only: a particle's own offset uses its own beta, which is why the
! concatenation stores every entry beta rather than one number.

do is = 1, ns
  sl(is) = is - 1
  ct_slice(is) = (is - 1) * prm%slice_spacing
enddo
t_slice = -ct_slice / c_light
z_slice = prm%beta0 * ct_slice

call H5Gcreate_f (f_id, 'coords', g_id, h5_err)
call group_note (f_id, 'coords', 'axis', 'Every axis of this file, once each.', err)

call fel_h5_int (g_id, 'record', '1', 'record', &
      'Record number, the axis every per-record dataset runs over.', 'record', rec, err)
call fel_h5_real (g_id, 's', 'm', 's', &
      'Path length along the line at each record, Bmad''s s. A VARIABLE on the record ' // &
      'axis: it repeats wherever two records land on one plane, so it cannot index.', &
      'record', stats%z(1:ir), err)
call fel_h5_int (g_id, 'ix_ele', '1', 'element', &
      'Lattice element holding each record. Its values index the ele axis, so it is ' // &
      'the join key onto the lattice/ table.', 'record', stats%ix_ele(1:ir), err)
call fel_h5_flag (g_id, 'at_element_end', 'at element end', &
      'One where the record is an element end. Selects the element_end axis, in order.', &
      'record', stats%at_end(1:ir), err)
call fel_h5_int (g_id, 'element_end', '1', 'element', &
      'The element index at each element end, the element_end axis itself.', &
      'element_end', e_ix, err)
call fel_h5_real (g_id, 's_element_end', 'm', 's at element end', &
      'Path length at each element end, which is coords/s where at_element_end is one.', &
      'element_end', e_s, err)
call fel_h5_int (g_id, 'slice', '1', 'slice', &
      'Slice number, the axis every per-slice dataset runs over.', 'slice', sl, err)
call fel_h5_real (g_id, 'ct_slice', 'm', 'ct', &
      'Light-travel distance of each slice ahead of the reference, window_sample ' // &
      'wavelengths per slice. EXACT and free of beta, and the coordinate slippage ' // &
      'counts in.', 'slice', ct_slice, err)
call fel_h5_real (g_id, 't_slice', 's', 't', &
      'Arrival time of each slice relative to the reference, -ct_slice/c exactly. More ' // &
      'negative toward the window head.', 'slice', t_slice, err)
call fel_h5_real (g_id, 'z_slice', 'm', 'z', &
      'Bmad z of each slice AT THE REFERENCE beta, beta0*ct_slice. A given particle''s ' // &
      'own offset uses its own beta, so this one number cannot serve for all of them.', &
      'slice', z_slice, err)
call fel_h5_str (g_id, 'bmad', 'phase space', &
      'Bmad phase-space coordinate names: the label axis of a centroid and the first ' // &
      'axis of a sigma matrix.', 'bmad', fel_stats_bmad$, err)
call fel_h5_str (g_id, 'bmad_col', 'phase space', &
      'The second axis of a Bmad sigma matrix, the same coordinates as bmad. Named ' // &
      'apart from it so that selecting one entry of a square matrix needs no rule.', &
      'bmad_col', fel_stats_bmad$, err)
call fel_h5_str (g_id, 'wavefront', 'wavefront', &
      'Wavefront moment coordinate names: the label axis of a field centroid and the ' // &
      'first axis of a field sigma.', 'wavefront', fel_stats_wf$, err)
call fel_h5_str (g_id, 'wavefront_col', 'wavefront', &
      'The second axis of a field sigma, the same coordinates as wavefront.', &
      'wavefront_col', fel_stats_wf$, err)
call fel_h5_str (g_id, 'bmad_unit', 'unit', &
      'The unit of each bmad coordinate. A dataset over that axis names it through ' // &
      '@unit_of_axis rather than spelling the list out where nothing can parse it.', &
      'bmad', fel_stats_bmad_u$, err)
call fel_h5_str (g_id, 'wavefront_unit', 'unit', &
      'The unit of each wavefront coordinate. See coords/bmad_unit.', &
      'wavefront', fel_stats_wf_u$, err)
call fel_h5_str (g_id, 'bmad_t', 'phase space + t', &
      'The six phase-space names plus t: the label axis of the per-slice extremes ' // &
      'rel_max and rel_min, which are seven-vectors in bunch_params_struct.', &
      'bmad_t', fel_stats_bmad_t$, err)
call fel_h5_str (g_id, 'bmad_t_unit', 'unit', &
      'The unit of each bmad_t entry. See coords/bmad_unit.', 'bmad_t', &
      fel_stats_bmad_t_u$, err)
call fel_h5_str (g_id, 'plane', 'plane', &
      'The three PROJECTED twiss planes, bunch_params_struct''s x, y and z.', 'plane', &
      fel_stats_planes$, err)
call fel_h5_str (g_id, 'mode', 'mode', &
      'The three NORMAL modes, bunch_params_struct''s a, b and c. A separate axis from ' // &
      'plane because the two are different decompositions of one beam, so nothing may ' // &
      'average across them. The labels are eigenvector-identified rather than ' // &
      'magnitude-sorted, which is why the harness compares mode emittances as a set.', &
      'mode', fel_stats_norm$, err)
if (err) return

! Which variable a viewer should put on the axis when it plots against the record or the
! element end. Both are indices, both have a path length riding on them as a variable,
! and only the writer knows which variable that is.

call H5LTset_attribute_string_f (g_id, 'record', 'plot_against', 's', h5_err)
err = err .or. (h5_err < 0)
call H5LTset_attribute_string_f (g_id, 'element_end', 'plot_against', 's_element_end', h5_err)
err = err .or. (h5_err < 0)
call H5LTset_attribute_string_f (g_id, 'slice', 'plot_against', 't_slice', h5_err)
err = err .or. (h5_err < 0)

! The join key and the mask say MACHINE-READABLY what their descriptions say in prose
! (bmad-stats R18, R19): ix_ele's values index the ele axis, and at_element_end selects
! the element_end axis's entries, in order.

call H5LTset_attribute_string_f (g_id, 'ix_ele', 'indexes', 'ele', h5_err)
err = err .or. (h5_err < 0)
call H5LTset_attribute_string_f (g_id, 'at_element_end', 'selects', 'element_end', h5_err)
err = err .or. (h5_err < 0)
if (err) return

! The head convention, which no per-slice plot can be drawn without: the migration
! invariant z_global = z_local + beta*(islice-1)*spacing puts the HEAD at the high slice
! index. The attribute states where the head sits on the INDEX, so both coordinates of
! that axis carry the same value while each description gives its own direction. A
! convention stated on one of a pair invites the reader to assume the pair agrees.

call H5LTset_attribute_string_f (g_id, 'ct_slice', 'head_direction', '+index', h5_err)
err = err .or. (h5_err < 0)
call H5LTset_attribute_string_f (g_id, 't_slice', 'head_direction', '+index', h5_err)
err = err .or. (h5_err < 0)
call H5LTset_attribute_string_f (g_id, 'z_slice', 'head_direction', '+index', h5_err)
err = err .or. (h5_err < 0)
call H5Gclose_f (g_id, h5_err)
if (err) return

! ------------------------------------------------------------------
! run/: the scalars the RUN produced, apart from params/, which holds what the user
! set and is written by fel_io_mod (one subgroup per honored input struct). The three
! counts restate axis lengths on purpose: n_element_end comes from the ACCUMULATOR's
! counter, so the harness checking it against the mask tests the walk's bookkeeping,
! which coords/element_end, packed from that mask, cannot.

call H5Gcreate_f (f_id, 'run', g_id, h5_err)
call group_note (f_id, 'run', 'run', 'What the run produced, one scalar each.', err)
call fel_h5_real (g_id, 'p0c', 'eV', 'p0c', 'Reference momentum times c.', '', stats%p0c, err)
call fel_h5_str (g_id, 'species', 'species', 'Particle species.', '', [prm%species], err)
call fel_h5_real (g_id, 'slice_spacing', 'm', 'slice spacing', &
      'Slice spacing, window_sample * lambda0. A LIGHT-TRAVEL distance, c times the ' // &
      'slice time separation: the grid is exactly uniform in t and ct, which is what ' // &
      'makes slippage a whole-slice shift. In Bmad z the separation is beta*this, per ' // &
      'particle. See coords/ct_slice.', '', prm%slice_spacing, err)
call fel_h5_int (g_id, 'n_record', '1', 'records', &
      'Records taken, the length of the record axis.', '', ir, err)
call fel_h5_int (g_id, 'n_element_end', '1', 'element ends', &
      'Element ends taken, the length of the element_end axis. From the accumulator''s ' // &
      'counter, and equal to the count of coords/at_element_end.', '', ie, err)
call fel_h5_int (g_id, 'n_slice', '1', 'slices', &
      'Slices in the time window, the length of the slice axis.', '', ns, err)
call H5Gclose_f (g_id, h5_err)
if (err) return

! ------------------------------------------------------------------
! beam/slice/: the per-record sufficient statistics, bunch_params_struct names, plus
! the derived quantities every consumer would otherwise re-derive.

allocate (cur(ns, ir), energy(ns, ir), sig_energy(ns, ir))
cur = c_light * stats%charge_live(:,1:ir) / prm%slice_spacing
energy = sqrt((stats%p0c * (1 + stats%b_centroid(6,:,1:ir)))**2 + m_electron**2)
! sigma_E = beta * p0c * sigma_pz, since dE/dp = p/E = beta. The factor is 1 to a few
! parts in 1e9 at an X-ray FEL's energy and it is not 1 for everything this tracker
! could be pointed at, so it is here rather than dropped.
sig_energy = stats%p0c * sqrt(max(0.0_rp, stats%b_sigma(36,:,1:ir))) * &
             stats%p0c * (1 + stats%b_centroid(6,:,1:ir)) / energy

call H5Gcreate_f (f_id, 'beam', b_id, h5_err)
call group_note (f_id, 'beam', 'beam', 'The electron beam, per slice and projected.', err)

call H5Gcreate_f (b_id, 'slice', g_id, h5_err)
call group_note (b_id, 'slice', 'per_slice', &
      'Sufficient statistics per slice, at every record.', err)
call fel_h5_real (g_id, 'centroid', fel_stats_bmad_unit$, 'centroid', &
      'Weighted centroid, Bmad phase space.', 'record,slice,bmad', &
      stats%b_centroid(:,:,1:ir), err)
call unit_axis_note (g_id, 'centroid', 'bmad', 1, err)
call fel_h5_real (g_id, 'sigma', fel_stats_bmad_unit$ // ' squared', 'sigma', &
      'Second moments about the centroid, the 6x6 at its natural rank.', &
      'record,slice,bmad,bmad_col', reshape(stats%b_sigma(:,:,1:ir), [6, 6, ns, ir]), err)
call unit_axis_note (g_id, 'sigma', 'bmad', 2, err)
call fel_h5_real (g_id, 'charge_live', 'C', 'charge', 'Live charge of the slice.', &
      'record,slice', stats%charge_live(:,1:ir), err)
call fel_h5_int (g_id, 'n_particle_live', '1', 'particles', &
      'Live macroparticles in the slice.', 'record,slice', stats%n_particle_live(:,1:ir), err)
call fel_h5_real (g_id, 't', 's', 't', &
      'Mean arrival time of the slice at this plane, relative to the reference.', &
      'record,slice', stats%t(:,1:ir), err)
call fel_h5_real (g_id, 'sigma_t', 's', 'sigma_t', &
      'Rms arrival-time spread within the slice.', 'record,slice', stats%sigma_t(:,1:ir), err)
call fel_h5_real (g_id, 'bunching', '1', 'bunching', &
      'Bunching |b| at the fundamental, charge weighted.', 'record,slice', &
      stats%bunching(:,1:ir), err)
call fel_h5_real (g_id, 'bunching_phase', 'rad', 'bunching phase', &
      'Bunching phase arg(b), carrying the run''s reference phase.', &
      'record,slice', stats%bunching_phase(:,1:ir), err)
call fel_h5_real (g_id, 'rel_max', 'm,1,m,1,m,1,s', 'max - centroid', &
      'Per-coordinate maximum over the slice RELATIVE TO THE CENTROID, ' // &
      'bunch_params_struct''s rel_max: order statistics no moment can reconstruct. ' // &
      'The envelope is centroid + rel_max. NaN for an empty slice.', &
      'record,slice,bmad_t', stats%b_rel_max(:,:,1:ir), err)
call unit_axis_note (g_id, 'rel_max', 'bmad_t', 1, err)
call fel_h5_real (g_id, 'rel_min', 'm,1,m,1,m,1,s', 'min - centroid', &
      'Per-coordinate minimum over the slice relative to the centroid. See rel_max.', &
      'record,slice,bmad_t', stats%b_rel_min(:,:,1:ir), err)
call unit_axis_note (g_id, 'rel_min', 'bmad_t', 1, err)

! The conveniences every consumer would otherwise re-derive, marked as the pure
! functions they are (bmad-stats R38). Ambient scalars (slice_spacing, p0c) are not
! named: the formula is in each description.

call fel_h5_real (g_id, 'current', 'A', 'current', &
      'Slice current, c * charge_live / slice_spacing.', 'record,slice', cur, err)
call fel_h5_dset_attr_strs (g_id, 'current', 'derived_from', [character(11):: 'charge_live'], err)
call fel_h5_real (g_id, 'energy', 'eV', 'energy', &
      'Mean total energy of the slice, Bmad''s convention (energy is eV, never gamma).', &
      'record,slice', energy, err)
call fel_h5_dset_attr_strs (g_id, 'energy', 'derived_from', [character(8):: 'centroid'], err)
call fel_h5_real (g_id, 'sigma_energy', 'eV', 'sigma_E', &
      'Rms energy spread of the slice.', 'record,slice', sig_energy, err)
call fel_h5_dset_attr_strs (g_id, 'sigma_energy', 'derived_from', &
      [character(8):: 'sigma', 'centroid'], err)
call H5Gclose_f (g_id, h5_err)
if (err) return

! beam/slice_twiss/ and beam/bunch/: the evaluated bunch_params, on the ELEMENT-END
! axis. The twiss planes are an AXIS, not six groups: a group named z cannot sit beside
! a z coordinate, which xarray and netCDF both refuse, and one array per quantity is 18
! datasets where six groups per set were 108. TWO axes, though, not one: twiss/ over the
! three projected planes and modes/ over the three normal modes, since averaging a
! projected emittance against an eigen-emittance is meaningless. They are subgroups
! rather than datasets here because one of the nine is named sigma, as
! bunch_params_struct names it, and the covariance matrix beside them is sigma too.

call H5Gcreate_f (b_id, 'slice_twiss', g_id, h5_err)
call group_note (b_id, 'slice_twiss', 'per_slice', &
      'Bmad''s evaluated bunch_params per slice, at element ends.', err)
call write_bp_slices (g_id, stats%e_slice(:, :, 1:ie), err)
call H5Gclose_f (g_id, h5_err)
if (err) return

call H5Gcreate_f (b_id, 'bunch', g_id, h5_err)
call group_note (b_id, 'bunch', 'projected', &
      'Bmad''s evaluated bunch_params for the whole window, at element ends. Pooled ' // &
      'from beam/slice''s moments by the covariance identity, so it is derived.', err)
call string_note (g_id, 'derived_from', [character(5):: 'slice'], err)
call write_bp_bunch (g_id, stats%e_bunch(:, 1:ie), err)
call H5Gclose_f (g_id, h5_err)
call H5Gclose_f (b_id, h5_err)
if (err) return

! ------------------------------------------------------------------
! field/: total/ always, then one group per component and per harmonic, all carrying
! the same dataset names. The group attributes say WHICH children are components and
! which are derived, so a reader summing children cannot double-count.

call H5Gcreate_f (f_id, 'field', b_id, h5_err)
call group_note (f_id, 'field', 'field', &
      'The radiation, per wavelength and per polarization component.', err)

if (allocated(stats%f2_power)) then
  call string_note (b_id, 'components', fel_stats_pol$, err)
  call write_field_total (b_id, stats%f_power(:,1:ir) + stats%f2_power(:,1:ir), &
                          stats%f_energy(:,1:ir) + stats%f2_energy(:,1:ir), &
                          stats%f_on_axis(:,1:ir) + stats%f2_on_axis(:,1:ir), 1, 2, err)
else
  call string_note (b_id, 'components', fel_stats_pol$(1:1), err)
  call write_field_total (b_id, stats%f_power(:,1:ir), stats%f_energy(:,1:ir), &
                          stats%f_on_axis(:,1:ir), 1, 1, err)
endif
if (err) return

call H5Gcreate_f (b_id, 'x', g_id, h5_err)
call group_note (b_id, 'x', 'component', 'The x polarization component.', err)
call write_field_component (g_id, stats%f_centroid(:,:,1:ir), stats%f_sigma(:,:,1:ir), &
        stats%f_power(:,1:ir), stats%f_energy(:,1:ir), stats%f_on_axis(:,1:ir), &
        stats%f_emit_x(:,1:ir), stats%f_emit_y(:,1:ir), stats%f_angles_valid(:,1:ir), err)
call H5Gclose_f (g_id, h5_err)
if (err) return

if (allocated(stats%f2_power)) then
  call H5Gcreate_f (b_id, 'y', g_id, h5_err)
  call group_note (b_id, 'y', 'component', 'The y polarization component.', err)
  call write_field_component (g_id, stats%f2_centroid(:,:,1:ir), stats%f2_sigma(:,:,1:ir), &
          stats%f2_power(:,1:ir), stats%f2_energy(:,1:ir), stats%f2_on_axis(:,1:ir), &
          stats%f2_emit_x(:,1:ir), stats%f2_emit_y(:,1:ir), stats%f_angles_valid(:,1:ir), err)
  call H5Gclose_f (g_id, h5_err)
  if (err) return
endif

! Harmonics: a group each, with its own total, since a detector separates colors and
! nothing may sum across wavelengths.

if (allocated(stats%fh_power)) then
  call H5LTset_attribute_int_f (f_id, 'field', 'harmonics', stats%fh_harm, &
                                int(size(stats%fh_harm), size_t), h5_err)
  err = err .or. (h5_err < 0)
  do ihh = 1, size(stats%fh_harm)
    write (hname, '(a, i0)') 'harm', stats%fh_harm(ihh)
    call H5Gcreate_f (b_id, trim(hname), s_id, h5_err)
    call group_note (b_id, trim(hname), 'harmonic', 'One harmonic field of the set.', err)
    call string_note (s_id, 'components', fel_stats_pol$(1:1), err)
    call fel_h5_attr_int (s_id, 'harmonic', stats%fh_harm(ihh), err)
    call write_field_total (s_id, stats%fh_power(:,1:ir,ihh), stats%fh_energy(:,1:ir,ihh), &
                            stats%fh_on_axis(:,1:ir,ihh), stats%fh_harm(ihh), 1, err)
    call H5Gcreate_f (s_id, 'x', g_id, h5_err)
    call group_note (s_id, 'x', 'component', 'The x polarization component.', err)
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
! Subroutine group_note (parent, name, kind, descrip, err)
!
! Routine to say what a group is: @kind for a reader that sorts them, @description for
! a reader who reads. Groups carried neither before, so beam/slice against beam/bunch
! had to be inferred from the name.
!-

subroutine group_note (parent, name, kind, descrip, err)

integer(hid_t) parent
integer h5e
logical err
character(*) name, kind, descrip

!

call H5LTset_attribute_string_f (parent, name, 'kind', kind, h5e)
err = err .or. (h5e < 0)
call H5LTset_attribute_string_f (parent, name, 'description', descrip, h5e)
err = err .or. (h5e < 0)

end subroutine group_note

!+
! Subroutine string_note (parent, name, attrib, val, err)
!
! Routine to set one string attribute on a group, reached through its parent by name.
! @components names the polarizations a field group holds and @derived_from what a
! derived group sums, so the always-written total/ cannot be taken for a component and
! then summed twice.
!
!-

subroutine string_note (obj_id, attrib, val, err)

integer(hid_t) obj_id
logical err, merr
character(*) attrib, val(:)

!

! An ARRAY, length one included, because these hold a LIST. A one-component file then
! parses exactly like a two-component one, where a bare string and a comma list would
! read as different types to anything that did not know to split.

call hdf5_write_attribute_string_rank1 (obj_id, attrib, val, merr)
err = err .or. merr

end subroutine string_note

!+
! Subroutine unit_axis_note (id, name, family, power, err)
!
! Routine to give a dataset whose entries have DIFFERENT units a machine-readable one
! beside the human string. That string is a comma list per coordinate, sometimes with
! "squared" on the end, which nothing can parse, so @unit_of_axis names the axis whose
! coords/<axis>_unit carries the units and @unit_power says how many factors of them one
! entry holds: 1 for a centroid, 2 for a second moment. @unit_power is a true scalar,
! since it is one number: see fel_h5_attr_int.
!-

subroutine unit_axis_note (id, name, family, power, err)

integer(hid_t) id
integer h5e, power
logical err
character(*) name, family

!

call H5LTset_attribute_string_f (id, name, 'unit_of_axis', family, h5e)
err = err .or. (h5e < 0)
call fel_h5_dset_attr_int (id, name, 'unit_power', power, err)

end subroutine unit_axis_note

!+
! Subroutine write_field_total (id, pow, ene, onax, harm, comps, err)
!
! Routine to write the total/ group of one wavelength: the sum over the live
! polarizations. ALWAYS written, whether one polarization is live or two, so that no
! reader has to ask what else the file holds before it knows what power means, and
! marked @derived_from so that summing the children of field/ cannot double-count.
!-

subroutine write_field_total (id, pow, ene, onax, harm, ncomp, err)

integer(hid_t) id, tt_id
integer harm, ncomp, h5e
real(rp) pow(:,:), ene(:,:), onax(:,:)
logical err
character(70) note

!

write (note, '(a, i0, a)') ' Summed over the live polarizations of harmonic ', harm, '.'

call H5Gcreate_f (id, 'total', tt_id, h5e)
call group_note (id, 'total', 'derived', 'The sum over this wavelength''s components.', err)
call string_note (tt_id, 'derived_from', fel_stats_pol$(1:ncomp), err)

call fel_h5_real (tt_id, 'power', 'W', 'power', 'Radiation power of the slice.' // trim(note), &
      'record,slice', pow, err)
call fel_h5_real (tt_id, 'energy', 'J', 'field energy', &
      'Field energy of the slice.' // trim(note), 'record,slice', ene, err)
call fel_h5_real (tt_id, 'on_axis_intensity', 'W/m^2', 'on-axis intensity', &
      'Intensity at the grid center.' // trim(note), 'record,slice', onax, err)
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

call fel_h5_real (id, 'centroid', fel_stats_wf_unit$, 'centroid', &
      'Intensity-weighted (x, theta_x, y, theta_y). The theta entries are NaN where ' // &
      'the angle moments were not computed.', 'record,slice,wavefront', cen, err)
call unit_axis_note (id, 'centroid', 'wavefront', 1, err)
call fel_h5_real (id, 'sigma', fel_stats_wf_unit$ // ' squared', 'sigma', &
      'Wigner second moments, the 4x4 at its natural rank. The theta rows are NaN ' // &
      'where the angle moments were not computed.', &
      'record,slice,wavefront,wavefront_col', &
      reshape(sig, [4, 4, size(sig,2), size(sig,3)]), err)
call unit_axis_note (id, 'sigma', 'wavefront', 2, err)
call fel_h5_real (id, 'power', 'W', 'power', 'Radiation power of this component.', &
      'record,slice', pow, err)
call fel_h5_real (id, 'energy', 'J', 'field energy', 'Field energy of this component.', &
      'record,slice', ene, err)
call fel_h5_real (id, 'on_axis_intensity', 'W/m^2', 'on-axis intensity', &
      'Intensity of this component at the grid center.', 'record,slice', onax, err)
call fel_h5_real (id, 'emit_x', 'm rad', 'emit_x', &
      'sqrt(det sigma_x-plane) = M^2 lambda / 4 pi. NaN without angle moments.', &
      'record,slice', ex, err)
call fel_h5_real (id, 'emit_y', 'm rad', 'emit_y', &
      'sqrt(det sigma_y-plane) = M^2 lambda / 4 pi. NaN without angle moments.', &
      'record,slice', ey, err)
call fel_h5_dset_attr_strs (id, 'emit_x', 'derived_from', [character(5):: 'sigma'], err)
call fel_h5_dset_attr_strs (id, 'emit_y', 'derived_from', [character(5):: 'sigma'], err)
call fel_h5_flag (id, 'angle_moments_valid', 'angles valid', &
      'One where the theta moments were computed (they cost three FFTs, so element ' // &
      'ends only).', 'record,slice', valid, err)

end subroutine write_field_component

!+
! Subroutine write_bp_bunch (id, rows, err)
!
! Routine to write the whole-window bunch_params rows: the packed 99-number row
! unpacked into named datasets, the twiss quantities over the plane axis.
!-

subroutine write_bp_bunch (id, rows, err)

integer(hid_t) id
real(rp) rows(:,:)
integer nr
logical err

!

nr = size(rows, 2)
call fel_h5_real (id, 'centroid', fel_stats_bmad_unit$, 'centroid', &
      'Weighted centroid of the whole window, Bmad phase space.', &
      'element_end,bmad', rows(1:6, :), err)
call unit_axis_note (id, 'centroid', 'bmad', 1, err)
call fel_h5_real (id, 'sigma', fel_stats_bmad_unit$ // ' squared', 'sigma', &
      'Second moments of the whole window, the 6x6 at its natural rank.', &
      'element_end,bmad,bmad_col', reshape(rows(7:42, :), [6, 6, nr]), err)
call unit_axis_note (id, 'sigma', 'bmad', 2, err)
call fel_h5_real (id, 'charge_live', 'C', 'charge', 'Live charge of the whole window.', &
      'element_end', rows(43, :), err)
call fel_h5_int (id, 'n_particle_live', '1', 'particles', &
      'Live macroparticles in the whole window.', 'element_end', nint(rows(44, :)), err)
call fel_h5_flag (id, 'twiss_valid', 'twiss valid', &
      'One where Bmad evaluated the twiss (it needs six live particles).', &
      'element_end', nint(rows(45, :)), err)
if (err) return

! Packed row layout: after the 45 leading numbers come six blocks of nine, one per
! decomposition entry in fel_stats_modes$ order, so one quantity is a stride-9 gather.
! Blocks 1 to 3 are the projected planes and 4 to 6 the normal modes, and they go to
! SEPARATE groups over separate axes: an eigen-emittance and a projected emittance are
! different quantities, and one axis holding both invites an average across them.

call write_twiss_group (id, 'twiss', 'plane', 0, rows, err)
call write_twiss_group (id, 'modes', 'mode', 3, rows, err)

end subroutine write_bp_bunch

!+
! Subroutine write_twiss_group (id, gname, axis, off, rows, err)
!
! Routine to write one decomposition's nine twiss quantities as nine datasets over its
! own axis, gathered from the packed rows starting at block off+1.
!-

subroutine write_twiss_group (id, gname, axis, off, rows, err)

integer(hid_t) id, t_id
real(rp) rows(:,:)
real(rp), allocatable :: tw(:,:)
integer off, im, k, jp, nr, h5e
logical err
character(*) gname, axis

!

nr = size(rows, 2)
call H5Gcreate_f (id, gname, t_id, h5e)
call group_note (id, gname, gname, 'The nine twiss quantities over the ' // axis // &
      ' axis, for the whole window. Pure functions of the sibling moments.', err)
call string_note (t_id, 'derived_from', [character(8):: 'centroid', 'sigma'], err)

allocate (tw(3, nr))
do jp = 1, 9
  k = 45 + 9 * off
  do im = 1, 3
    tw(im, :) = rows(k + jp, :)
    k = k + 9
  enddo
  call fel_h5_real (t_id, trim(fel_stats_twiss$(jp)), trim(fel_stats_tunit$(jp)), &
        trim(fel_stats_twiss$(jp)), 'Twiss ' // trim(fel_stats_twiss$(jp)) // &
        ' of the whole window, per ' // axis // '. NaN where twiss_valid is zero.', &
        'element_end,' // axis, tw, err)
  if (err) exit
enddo
call H5Gclose_f (t_id, h5e)

end subroutine write_twiss_group

!+
! Subroutine write_bp_slices (id, rows, err)
!
! Routine to write the per-slice bunch_params rows, the same names as write_bp_bunch
! with the slice axis added.
!-

subroutine write_bp_slices (id, rows, err)

integer(hid_t) id
real(rp) rows(:,:,:)
integer nsl, nr
logical err

!

nsl = size(rows, 2);  nr = size(rows, 3)
call fel_h5_real (id, 'centroid', fel_stats_bmad_unit$, 'centroid', &
      'Weighted centroid per slice, Bmad phase space.', &
      'element_end,slice,bmad', rows(1:6, :, :), err)
call unit_axis_note (id, 'centroid', 'bmad', 1, err)
call fel_h5_real (id, 'sigma', fel_stats_bmad_unit$ // ' squared', 'sigma', &
      'Second moments per slice, the 6x6 at its natural rank.', &
      'element_end,slice,bmad,bmad_col', reshape(rows(7:42, :, :), [6, 6, nsl, nr]), err)
call unit_axis_note (id, 'sigma', 'bmad', 2, err)
call fel_h5_real (id, 'charge_live', 'C', 'charge', 'Live charge of the slice.', &
      'element_end,slice', rows(43, :, :), err)
call fel_h5_int (id, 'n_particle_live', '1', 'particles', &
      'Live macroparticles in the slice.', 'element_end,slice', nint(rows(44, :, :)), err)
call fel_h5_flag (id, 'twiss_valid', 'twiss valid', &
      'One where Bmad evaluated the twiss (it needs six live particles).', &
      'element_end,slice', nint(rows(45, :, :)), err)
if (err) return

! Two groups over two axes: see write_bp_bunch.

call write_slice_twiss_group (id, 'twiss', 'plane', 0, rows, err)
call write_slice_twiss_group (id, 'modes', 'mode', 3, rows, err)

end subroutine write_bp_slices

!+
! Subroutine write_slice_twiss_group (id, gname, axis, off, rows, err)
!
! Routine to write one decomposition's nine twiss quantities per slice, the same names
! as write_twiss_group with the slice axis added.
!-

subroutine write_slice_twiss_group (id, gname, axis, off, rows, err)

integer(hid_t) id, t_id
real(rp) rows(:,:,:)
real(rp), allocatable :: tw(:,:,:)
integer off, im, k, jp, nsl, nr, h5e
logical err
character(*) gname, axis

!

nsl = size(rows, 2);  nr = size(rows, 3)
call H5Gcreate_f (id, gname, t_id, h5e)
call group_note (id, gname, gname, 'The nine twiss quantities over the ' // axis // &
      ' axis, per slice. Pure functions of the sibling moments.', err)
call string_note (t_id, 'derived_from', [character(8):: 'centroid', 'sigma'], err)

allocate (tw(3, nsl, nr))
do jp = 1, 9
  k = 45 + 9 * off
  do im = 1, 3
    tw(im, :, :) = rows(k + jp, :, :)
    k = k + 9
  enddo
  call fel_h5_real (t_id, trim(fel_stats_twiss$(jp)), trim(fel_stats_tunit$(jp)), &
        trim(fel_stats_twiss$(jp)), 'Twiss ' // trim(fel_stats_twiss$(jp)) // &
        ' per slice and ' // axis // '. NaN where twiss_valid is zero.', &
        'element_end,slice,' // axis, tw, err)
  if (err) exit
enddo
call H5Gclose_f (t_id, h5e)

end subroutine write_slice_twiss_group

end subroutine fel_stats_write

end module fel_stats_mod
