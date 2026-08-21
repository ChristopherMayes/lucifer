!+
! Module fel_io_mod
!
! The file artifacts of a run: field-set dumps (Genesis and openPMD formats), the
! openPMD field imports, the escaped-slice bank and its drain, the full-pulse
! reconstruction, the stats-file finalize, and the wake-eloss block writer. Every
! procedure takes the run state (fel_run_struct) explicitly, and NOTHING HERE STOPS:
! errors return through err_flag and the caller decides -- the library contract
! (manual sec:program). The print lines are unchanged from when this code lived in
! the driver; the check scripts parse them.
!-

module fel_io_mod

use fel_struct
use wavefront_hdf5_mod
use wavefront_openpmd_mod

implicit none

contains

!------------------------------------------------------------------------------
! Write the whole field set at the given filename prefix, honoring wavefront_format:
! Genesis dumps (one polarization per file: -x/-y when Ey is live) and/or openPMD
! EXT_Wavefront (both polarizations as components of ONE mesh record). The
! fundamental keeps the pre-harmonic names; a harmonic's files carry -h<h>. Every
! record is unrotated to time order first (each field owns its rotation state; they
! move in lockstep, but the cshift must run per record).

subroutine fel_dump_field_set (run, prefix, err_flag)

type (fel_run_struct), target :: run
type (fel_field_struct), pointer :: ffield(:)
character(*) prefix
logical err_flag
integer ihh
logical eerr
character(8) hsuf
character(8) wavefront_format
real(rp), pointer :: z_now
integer n_harm

!

err_flag = .true.
ffield => run%ffield
z_now => run%z_now
wavefront_format = run%global%wavefront_format
n_harm = run%n_harm

do ihh = 1, n_harm
  if (ffield(ihh)%slip%first /= 0) then
    ffield(ihh)%wf%Ex = cshift(ffield(ihh)%wf%Ex, shift = ffield(ihh)%slip%first, dim = 3)
    if (allocated(ffield(ihh)%wf%Ey)) &
        ffield(ihh)%wf%Ey = cshift(ffield(ihh)%wf%Ey, shift = ffield(ihh)%slip%first, dim = 3)
    ffield(ihh)%slip%first = 0
  endif

  hsuf = ''
  if (ffield(ihh)%harm /= 1) write (hsuf, '(a, i0)') '-h', ffield(ihh)%harm

  if (wavefront_format /= 'openpmd') then         ! genesis or both
    if (allocated(ffield(ihh)%wf%Ey)) then        ! One component per Genesis file.
      call wavefront_write_genesis4 (ffield(ihh)%wf, prefix // trim(hsuf) // '-x.fld.h5', eerr, 'x')
      if (eerr) stop 1
      call wavefront_write_genesis4 (ffield(ihh)%wf, prefix // trim(hsuf) // '-y.fld.h5', eerr, 'y')
      if (eerr) stop 1
      print '(a)', '  ' // prefix // trim(hsuf) // '-{x,y}.fld.h5'
    else
      call wavefront_write_genesis4 (ffield(ihh)%wf, prefix // trim(hsuf) // '.fld.h5', eerr, 'x')
      if (eerr) stop 1
      print '(a)', '  ' // prefix // trim(hsuf) // '.fld.h5'
    endif
  endif

  if (wavefront_format /= 'genesis') then         ! openpmd or both
    call wavefront_write_openpmd (ffield(ihh)%wf, prefix // trim(hsuf) // '.wf.h5', z_now, eerr)
    if (eerr) stop 1
    print '(a)', '  ' // prefix // trim(hsuf) // '.wf.h5'
  endif
enddo

err_flag = .false.

end subroutine fel_dump_field_set

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Read an openPMD wavefront into field-set entry ihh (the fundamental import path).
! The photon energy must be the fundamental's -- a file carrying a harmonic in
! field_file(1) is refused by name.

subroutine fel_read_openpmd_into_field (run, fname, ihh, err_flag)

type (fel_run_struct), target :: run
type (fel_field_struct), pointer :: ffield(:)
type (fel_beam_struct), pointer :: fbeam
character(*) fname
integer ihh
real(rp) e_photon
logical rerr, err_flag

!

err_flag = .true.
ffield => run%ffield
fbeam => run%fbeam

call wavefront_read_openpmd (ffield(ihh)%wf, fname, rerr, e_photon)
if (rerr) stop 1
if (abs(ffield(ihh)%wf%wavelength - fbeam%wavelength) > 1e-6_rp * fbeam%wavelength) then
  print '(a)', 'fel_track_test: the openPMD file in field_file(1) does not carry the FUNDAMENTAL:'
  print '(a, 2es20.12)', '  its photonEnergy wavelength vs the beam: ', &
                         ffield(ihh)%wf%wavelength, fbeam%wavelength
  err_flag = .true.;  return
endif
ffield(ihh)%wf%dz = fbeam%slice_spacing
ffield(ihh)%wf%wavelength = fbeam%wavelength   ! One wavelength authority (the 1e-12
                                               ! beam/field consistency check's spirit).

err_flag = .false.

end subroutine fel_read_openpmd_into_field

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Import a harmonic field: match the file's photon energy to the field-set entry
! carrying that harmonic (no match is refused by name), require the fundamental's
! grid and window, and keep the walk's wavelength convention (fundamental / h).

subroutine fel_import_harmonic_field (run, fname, err_flag)

type (fel_run_struct), target :: run
type (fel_field_struct), pointer :: ffield(:)
type (fel_beam_struct), pointer :: fbeam
type (wavefront_struct), pointer :: wf
character(*) fname
type (wavefront_struct) wtmp
real(rp) e_photon, e1
integer ihh, n_harm
logical rerr, err_flag

!

err_flag = .true.
ffield => run%ffield
fbeam => run%fbeam
wf => run%ffield(1)%wf
n_harm = run%n_harm

call wavefront_read_openpmd (wtmp, fname, rerr, e_photon)
if (rerr) stop 1
e1 = h_planck * c_light / fbeam%wavelength * e_charge

do ihh = 2, n_harm
  if (abs(e_photon - ffield(ihh)%harm * e1) <= 1e-6_rp * ffield(ihh)%harm * e1) exit
enddo
if (ihh > n_harm) then
  print '(a)', 'fel_track_test: the photonEnergy of ' // trim(fname)
  print '(a, es13.5, a)', '  (', e_photon, ' J) matches NO field of this run''s harmonics list.'
  err_flag = .true.;  return
endif

if (size(wtmp%Ex,1) /= size(wf%Ex,1) .or. size(wtmp%Ex,3) /= size(wf%Ex,3) .or. &
    abs(wtmp%dx - wf%dx) > 1e-12_rp * wf%dx) then
  print '(a)', 'fel_track_test: harmonic import ' // trim(fname)
  print '(a)', '  does not match the fundamental''s grid and window (one time window, one grid).'
  err_flag = .true.;  return
endif

call move_alloc (wtmp%Ex, ffield(ihh)%wf%Ex)
if (allocated(wtmp%Ey)) call move_alloc (wtmp%Ey, ffield(ihh)%wf%Ey)
ffield(ihh)%wf%dx = wtmp%dx;  ffield(ihh)%wf%dy = wtmp%dy
ffield(ihh)%wf%dz = fbeam%slice_spacing
ffield(ihh)%wf%wavelength = wf%wavelength / ffield(ihh)%harm

err_flag = .false.

end subroutine fel_import_harmonic_field

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Drain the slippage bank: each transmitted slice streams to <out_root>-escaped.fld.h5
! as it leaves (peak memory a handful of grid planes), with its wavefront_params (full
! 4x4 at bank time -- exactly what the analytic free-space propagation of pulse
! statistics needs) and its transmission z. Genesis field-file conventions (dfl units)
! so existing tooling reads it; root datasets land at finalize.

subroutine fel_drain_bank (run, ihh, err_flag)

type (fel_run_struct), target :: run
type (fel_field_struct), pointer :: ffield(:)
type (fel_beam_struct), pointer :: fbeam
real(rp), pointer :: z_now
logical two_pol, err_flag
integer n_harm

type (wavefront_params_struct) pms
type (fel_bank_struct), pointer :: bnk
type (wavefront_struct), pointer :: wfl
integer(hid_t) g_id
real(rp), allocatable :: work(:), gz(:,:), gp(:,:,:)
real(rp) dfl_scale
integer ihh, k, nx, h5e
logical berr
character(20) gname

!

err_flag = .true.
ffield => run%ffield
fbeam => run%fbeam
z_now => run%z_now
two_pol = run%two_pol
n_harm = run%n_harm

bnk => ffield(ihh)%bank
wfl => ffield(ihh)%wf
if (bnk%n == 0) then
  err_flag = .false.
  return
endif
nx = size(wfl%Ex, 1)
dfl_scale = wfl%dx / sqrt(2 * (mu_0_vac * c_light))
allocate (work(nx*nx))

if (run%esc_id(ihh) == 0) then
  call hdf5_open_file (trim(fel_escaped_file_name(run, ihh)), 'WRITE', run%esc_id(ihh), berr)
  if (berr) stop 1
endif

do k = 1, bnk%n
  run%n_banked(ihh) = run%n_banked(ihh) + 1
  if (.not. allocated(run%bank_z)) then
    allocate (run%bank_z(1024, n_harm), run%bank_pms(25, 1024, n_harm))
  elseif (run%n_banked(ihh) > size(run%bank_z, 1)) then
    call move_alloc (run%bank_z, gz);  call move_alloc (run%bank_pms, gp)
    allocate (run%bank_z(2*size(gz,1), n_harm), run%bank_pms(25, 2*size(gz,1), n_harm))
    run%bank_z(1:size(gz,1), :) = gz;  run%bank_pms(:, 1:size(gz,1), :) = gp
    deallocate (gz, gp)
  endif

  call wavefront_params_of_plane (bnk%plane(:,:,k), wfl%dx, wfl%wavelength, &
                                  fbeam%slice_spacing, pms, .true., berr)
  if (berr) stop 1
  pms%s = z_now
  run%bank_z(run%n_banked(ihh), ihh) = z_now
  run%bank_pms(:, run%n_banked(ihh), ihh) = [pms%centroid, reshape(pms%sigma, [16]), &
                           pms%energy, pms%power, pms%on_axis_intensity, pms%emit_x, pms%emit_y]
  if (two_pol) then             ! The y component's params add to the banked energy.
    call wavefront_params_of_plane (bnk%plane_y(:,:,k), wfl%dx, wfl%wavelength, &
                                    fbeam%slice_spacing, pms, .true., berr)
    if (berr) stop 1
  endif

  write (gname, '(a, i0.6)') 'slice', run%n_banked(ihh)
  call H5Gcreate_f (run%esc_id(ihh), trim(gname), g_id, h5e)
  if (h5e < 0) stop 1
  work = dfl_scale * reshape(real(bnk%plane(:,:,k), rp), [nx*nx])
  call hdf5_write_dataset_real (g_id, 'field-real', work, berr);  if (berr) stop 1
  work = dfl_scale * reshape(aimag(bnk%plane(:,:,k)), [nx*nx])
  call hdf5_write_dataset_real (g_id, 'field-imag', work, berr);  if (berr) stop 1
  if (two_pol) then
    work = dfl_scale * reshape(real(bnk%plane_y(:,:,k), rp), [nx*nx])
    call hdf5_write_dataset_real (g_id, 'field-real-y', work, berr);  if (berr) stop 1
    work = dfl_scale * reshape(aimag(bnk%plane_y(:,:,k)), [nx*nx])
    call hdf5_write_dataset_real (g_id, 'field-imag-y', work, berr);  if (berr) stop 1
    call hdf5_write_dataset_real (g_id, 'wavefront_params_y', &
            [pms%centroid, reshape(pms%sigma, [16]), pms%energy, pms%power, &
             pms%on_axis_intensity, pms%emit_x, pms%emit_y], berr);  if (berr) stop 1
  endif
  call hdf5_write_dataset_real (g_id, 'z_transmit', [z_now], berr);  if (berr) stop 1
  call hdf5_write_dataset_real (g_id, 'wavefront_params', run%bank_pms(:, run%n_banked(ihh), ihh), berr);  if (berr) stop 1
  call H5Gclose_f (g_id, h5e)
enddo

err_flag = .false.

end subroutine fel_drain_bank

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Escaped-file name of field ihh: the fundamental keeps its pre-harmonic name.

function fel_escaped_file_name (run, ihh) result (fname)

type (fel_run_struct), target :: run
type (fel_field_struct), pointer :: ffield(:)
character(400) out_root
integer ihh
character(420) fname

ffield => run%ffield
out_root = run%global%out_root

if (ffield(ihh)%harm == 1) then
  fname = trim(out_root) // '-escaped.fld.h5'
else
  write (fname, '(2a, i0, a)') trim(out_root), '-escaped-h', ffield(ihh)%harm, '.fld.h5'
endif

end function fel_escaped_file_name

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Finalize: the stats file always; with keep_escaped_field also the escaped file's root
! datasets and the FULL PULSE at the exit plane -- each banked slice read back, free-
! space propagated over z_end - z_transmit (transmitted light is fixed information;
! undulator vacuum is free space for light with no beam under it), and concatenated
! above the live window: earliest-transmitted light is furthest ahead, so banked slice
! k lands at pulse index nslice + (run%n_banked - k + 1). The caller has already unrotated
! the live window (the final-dump cshift).

subroutine fel_finalize_diagnostics (run, err_flag)

type (fel_run_struct), target :: run
type (fel_field_struct), pointer :: ffield(:)
type (fel_beam_struct), pointer :: fbeam
integer nx, h5e, ihh, n_harm, nslice
logical ferr, err_flag, two_pol, keep_escaped_field
character(400) out_root
character(420) fname_h

!

err_flag = .true.
ffield => run%ffield
fbeam => run%fbeam
n_harm = run%n_harm
nslice = run%nslice
two_pol = run%two_pol
keep_escaped_field = run%global%keep_escaped_field
out_root = run%global%out_root

call fel_stats_write (run%stats, trim(out_root) // '.stats.h5', ferr)
if (ferr) stop 1
print '(a)', '  ' // trim(out_root) // '.stats.h5'

if (.not. keep_escaped_field) then
  err_flag = .false.
  return
endif

do ihh = 1, n_harm
  nx = size(ffield(ihh)%wf%Ex, 1)

  if (run%esc_id(ihh) /= 0) then
    call hdf5_write_dataset_int  (run%esc_id(ihh), 'gridpoints',   [nx],                        ferr)
    call hdf5_write_dataset_real (run%esc_id(ihh), 'gridsize',     [ffield(ihh)%wf%dx],         ferr)
    call hdf5_write_dataset_real (run%esc_id(ihh), 'refposition',  [ffield(ihh)%wf%ref_position], ferr)
    call hdf5_write_dataset_real (run%esc_id(ihh), 'wavelength',   [ffield(ihh)%wf%wavelength], ferr)
    call hdf5_write_dataset_int  (run%esc_id(ihh), 'slicecount',   [run%n_banked(ihh)],             ferr)
    call hdf5_write_dataset_real (run%esc_id(ihh), 'slicespacing', [fbeam%slice_spacing],       ferr)
    call H5Fclose_f (run%esc_id(ihh), h5e)
    run%esc_id(ihh) = 0
    print '(a)', '  ' // trim(fel_escaped_file_name(run, ihh))
  endif

  ! The full pulse at the exit plane; with two live polarizations, one file per
  ! component (Genesis's format holds one). Harmonic pulses carry -h<h>.

  if (two_pol) then
    call fel_write_pulse_file (run, trim(out_root) // '-pulse-x.fld.h5', .false., ihh, ferr)
    if (ferr) return
    call fel_write_pulse_file (run, trim(out_root) // '-pulse-y.fld.h5', .true., ihh, ferr)
    if (ferr) return
    print '(a)', '  ' // trim(out_root) // '-pulse-{x,y}.fld.h5'
  elseif (ffield(ihh)%harm == 1) then
    call fel_write_pulse_file (run, trim(out_root) // '-pulse.fld.h5', .false., ihh, ferr)
    if (ferr) return
    print '(a)', '  ' // trim(out_root) // '-pulse.fld.h5'
  else
    write (fname_h, '(2a, i0, a)') trim(out_root), '-pulse-h', ffield(ihh)%harm, '.fld.h5'
    call fel_write_pulse_file (run, trim(fname_h), .false., ihh, ferr)
    if (ferr) return
    print '(a)', '  ' // trim(fname_h)
  endif
enddo

err_flag = .false.

end subroutine fel_finalize_diagnostics

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! One component of the full exit-plane pulse (Genesis's format holds one per file):
! the live window's slices, then each banked slice read back from the escaped file
! and free-space propagated over z_end - z_transmit. Shared by the -x and -y files.

subroutine fel_write_pulse_file (run, fname, use_y, ihh, err_flag)

type (fel_run_struct), target :: run
type (fel_field_struct), pointer :: ffield(:)
type (fel_beam_struct), pointer :: fbeam
type (wavefront_struct) wf1
real(rp), pointer :: z_now
character(*) fname
logical use_y, err_flag
integer ihh, n_harm, nslice
type (wavefront_struct), pointer :: wfl
integer(hid_t) p_id, e_id, g_id
real(rp), allocatable :: work(:), re_w(:), im_w(:)
real(rp) dfl_scale
integer k, is_f, nx, h5e
logical ferr
character(24) gname, dset_r, dset_i

!

err_flag = .true.
ffield => run%ffield
fbeam => run%fbeam
z_now => run%z_now
n_harm = run%n_harm
nslice = run%nslice

wfl => ffield(ihh)%wf
nx = size(wfl%Ex, 1)
dfl_scale = wfl%dx / sqrt(2 * (mu_0_vac * c_light))
dset_r = 'field-real';  dset_i = 'field-imag'
if (use_y) then
  dset_r = 'field-real-y';  dset_i = 'field-imag-y'
endif

call hdf5_open_file (trim(fname), 'WRITE', p_id, ferr);  if (ferr) stop 1
call hdf5_write_dataset_int  (p_id, 'gridpoints',   [nx],                       ferr)
call hdf5_write_dataset_real (p_id, 'gridsize',     [wfl%dx],                   ferr)
call hdf5_write_dataset_real (p_id, 'refposition',  [wfl%ref_position],         ferr)
call hdf5_write_dataset_real (p_id, 'wavelength',   [wfl%wavelength],           ferr)
call hdf5_write_dataset_int  (p_id, 'slicecount',   [nslice + run%n_banked(ihh)],   ferr)
call hdf5_write_dataset_real (p_id, 'slicespacing', [fbeam%slice_spacing],      ferr)

allocate (work(nx*nx), re_w(nx*nx), im_w(nx*nx))

do is_f = 1, nslice
  write (gname, '(a, i0.6)') 'slice', is_f
  call H5Gcreate_f (p_id, trim(gname), g_id, h5e)
  if (use_y) then
    work = dfl_scale * reshape(real(wfl%Ey(:,:,is_f), rp), [nx*nx])
  else
    work = dfl_scale * reshape(real(wfl%Ex(:,:,is_f), rp), [nx*nx])
  endif
  call hdf5_write_dataset_real (g_id, 'field-real', work, ferr);  if (ferr) stop 1
  if (use_y) then
    work = dfl_scale * reshape(aimag(wfl%Ey(:,:,is_f)), [nx*nx])
  else
    work = dfl_scale * reshape(aimag(wfl%Ex(:,:,is_f)), [nx*nx])
  endif
  call hdf5_write_dataset_real (g_id, 'field-imag', work, ferr);  if (ferr) stop 1
  call H5Gclose_f (g_id, h5e)
enddo

if (run%n_banked(ihh) > 0) then
  if (.not. allocated(wf1%Ex)) then
    allocate (wf1%Ex(nx, nx, 1))
    wf1%dx = wfl%dx;  wf1%dy = wfl%dy;  wf1%dz = fbeam%slice_spacing
  endif
  wf1%wavelength = wfl%wavelength      ! Per field: banked light drifts at ITS wavelength.

  call hdf5_open_file (trim(fel_escaped_file_name(run, ihh)), 'READ', e_id, ferr);  if (ferr) stop 1
  do k = 1, run%n_banked(ihh)
    write (gname, '(a, i0.6)') 'slice', k
    g_id = hdf5_open_group (e_id, trim(gname), ferr, .true.);  if (ferr) stop 1
    call hdf5_read_dataset_real (g_id, trim(dset_r), re_w, ferr, trim(gname));  if (ferr) stop 1
    call hdf5_read_dataset_real (g_id, trim(dset_i), im_w, ferr, trim(gname));  if (ferr) stop 1
    call H5Gclose_f (g_id, h5e)
    wf1%Ex(:,:,1) = reshape(cmplx(re_w, im_w, wf_rp), [nx, nx]) / dfl_scale

    call wavefront_drift (wf1, z_now - run%bank_z(k, ihh), ferr);  if (ferr) stop 1

    write (gname, '(a, i0.6)') 'slice', nslice + (run%n_banked(ihh) - k + 1)
    call H5Gcreate_f (p_id, trim(gname), g_id, h5e)
    work = dfl_scale * reshape(real(wf1%Ex(:,:,1), rp), [nx*nx])
    call hdf5_write_dataset_real (g_id, 'field-real', work, ferr);  if (ferr) stop 1
    work = dfl_scale * reshape(aimag(wf1%Ex(:,:,1)), [nx*nx])
    call hdf5_write_dataset_real (g_id, 'field-imag', work, ferr);  if (ferr) stop 1
    call H5Gclose_f (g_id, h5e)
  enddo
  call H5Fclose_f (e_id, h5e)
endif

call H5Fclose_f (p_id, h5e)

err_flag = .false.

end subroutine fel_write_pulse_file

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine fel_write_wake_block (run, z)

type (fel_run_struct), target :: run

! One block of per-slice eloss, z-stamped. Written at the hoisted update and at every
! migration-stride recompute; the energy-bookkeeping and stale-wake checks parse these.

real(rp) z
integer is_w

write (run%iu_wake, '(a, es22.14)') '# z = ', z
do is_w = 1, run%nslice
  write (run%iu_wake, '(i8, es24.16)') is_w, run%coll%wake%eloss(is_w)
enddo

end subroutine fel_write_wake_block

end module fel_io_mod
