!+
! Module fel_io_mod
!
! Routines to write the file artifacts of a run: the field-set dumps (Genesis and openPMD
! formats), the openPMD field imports, the escaped-slice bank and its drain, the
! full-pulse reconstruction, the stats-file finalize, and the wake-eloss block writer.
!
! Every procedure takes the run state (fel_run_struct) explicitly, and nothing here
! stops: errors return through err_flag and the caller decides -- the library contract
! (manual sec:program). All terminal output goes through out_io. The file-name echo
! lines are s_blank$ so scripts see bare names.
!-

module fel_io_mod

use fel_struct
use fel_input_mod
use wavefront_openpmd_mod
!$ use omp_lib

implicit none

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_dump_beam (run, ele, prefix, err_flag)
!
! Routine to write the beam at the given filename prefix as openPMD (.beam.h5), the
! only particle format this tracker writes. Genesis .par.h5 conversion lives in
! lucifer/tests/scripts/convert_genesis.py, outside the physics.
!
! Input:
!   run       -- fel_run_struct: Run state.
!   ele       -- ele_struct: Element the beam sits at.
!   prefix    -- character(*): Filename prefix. Format suffixes are appended.
!
! Output:
!   err_flag  -- logical: Set True if a file could not be written. False otherwise.
!-

subroutine fel_dump_beam (run, ele, prefix, err_flag)

type (fel_run_struct), target :: run
type (ele_struct) ele
character(*) prefix
logical err_flag
logical eerr

!

err_flag = .true.

call fel_write_openpmd_beam (run%fbeam, ele, prefix // '.beam.h5', eerr)
if (eerr) return

err_flag = .false.

end subroutine fel_dump_beam

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_dump_field_set (run, prefix, err_flag)
!
! Routine to write the whole field set at the given filename prefix as openPMD
! EXT_Wavefront (.wf.h5), both polarizations as components of ONE mesh record. Genesis
! .fld.h5 conversion lives in lucifer/tests/scripts/convert_genesis.py.
! The fundamental keeps the pre-harmonic names. A harmonic's files carry -h<h>.
! Every record is unrotated to time order first: each field owns its rotation state.
! The fields move in lockstep, but the cshift must run per record.
!
! Input:
!   run       -- fel_run_struct: Run state with the field set to dump.
!   prefix    -- character(*): Filename prefix. Format-specific suffixes are appended.
!
! Output:
!   run       -- fel_run_struct: Field records unrotated to time order (slip%first = 0).
!   err_flag  -- logical: Set True if a file could not be written. False otherwise.
!-

subroutine fel_dump_field_set (run, prefix, err_flag)

type (fel_run_struct), target :: run
type (fel_field_struct), pointer :: ffield(:)
character(*) prefix
logical err_flag
integer ihh
logical eerr
character(8) hsuf
real(rp), pointer :: z_now
integer n_harm
character(*), parameter :: r_name = 'fel_dump_field_set'

!

err_flag = .true.
ffield => run%ffield
z_now => run%z_now
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

  call wavefront_write_openpmd (ffield(ihh)%wf, prefix // trim(hsuf) // '.wf.h5', z_now, eerr)
  if (eerr) return
enddo

err_flag = .false.

end subroutine fel_dump_field_set

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_write_header (run)
!
! Routine to print the run's configuration as one framed block for a human to read:
! what lattice, what beam, what radiation, what is switched on, where the output goes.
! This replaces a scatter of loose informational messages. Display only -- the same
! values, resolved, are echoed machine-readably into the stats file's Meta group.
!
! Input:
!   run -- fel_run_struct: The fully set-up run state.
!-

subroutine fel_write_header (run)

type (fel_run_struct), target :: run
type (branch_struct), pointer :: branch
character(200) line
character(*), parameter :: r_name = 'lucifer'
character(*), parameter :: rule = &
      '================================================================================'
integer nfel, n_omp

!

branch => run%lat%branch(0)
nfel = count(run%is_fel)
n_omp = 1
!$ n_omp = omp_get_max_threads()

call out_io (s_blank$, r_name, rule)
call out_io (s_blank$, r_name, ' Lucifer -- FEL tracking in Bmad')
call out_io (s_blank$, r_name, &
      '--------------------------------------------------------------------------------')

write (line, '(a, a)') ' Lattice     ', trim(run%lat_file)
call out_io (s_blank$, r_name, trim(line))
write (line, '(a, i0, a, f0.3, a, i0, a)') '             ', branch%n_ele_track, ' elements, ', &
      branch%ele(branch%n_ele_track)%s, ' m, ', nfel, ' FEL segments'
call out_io (s_blank$, r_name, trim(line))

if (run%nslice == 1) then
  write (line, '(a, i0, a, f0.2)') ' Beam        1 slice x ', run%fbeam%slice(1)%n, &
        ' particles, gamma0 = ', run%gamma0
else
  write (line, '(a, i0, a, i0, a, f0.2)') ' Beam        ', run%nslice, ' slices x ', &
        run%fbeam%slice(1)%n, ' particles per slice, gamma0 = ', run%gamma0
endif
call out_io (s_blank$, r_name, trim(line))

write (line, '(5a, i0, 2a)') ' Radiation   lambda0 = ', trim(adjustl(fel_si_str(run%winit%lambda0, 'm'))), &
      ', slice spacing ', trim(adjustl(fel_si_str(run%fbeam%slice_spacing, 'm'))), ', ', &
      run%n_harm, ' field(s), grid half width ', trim(adjustl(fel_si_str(run%winit%grid_half_width, 'm')))
call out_io (s_blank$, r_name, trim(line))

write (line, '(a, l1, a, l1, a, l1)') ' Switches    sr wakes ', run%coll%wake%on, &
      ', space charge ', run%coll%efield%on, ', radiation damping ', bmad_com%radiation_damping_on
call out_io (s_blank$, r_name, trim(line))

write (line, '(3a, i0)') ' Output      out_root "', trim(run%global%out_root), '", threads ', n_omp
call out_io (s_blank$, r_name, trim(line))
call out_io (s_blank$, r_name, rule)

end subroutine fel_write_header

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_write_footer (run)
!
! Routine to close a run for a human: how far it went and how long it took, what the
! light ended up as, and which files were written. The file names were formerly printed
! as the writers went. Collecting them here means the last thing on screen is the list
! of what to open.
!
! Input:
!   run -- fel_run_struct: The finished run state.
!-

subroutine fel_write_footer (run)

type (fel_run_struct), target :: run
type (fel_stats_struct), pointer :: stats
character(300) line
character(*), parameter :: r_name = 'lucifer'
integer ir, ih, n_listed
real(rp) pow, ene, bun
character(8) hsuf

!

stats => run%stats
call out_io (s_blank$, r_name, &
      '--------------------------------------------------------------------------------')

write (line, '(a, f0.3, a, i0, a)') ' Done        ', run%z_now, ' m, ', run%stats%iend, ' element ends'
call out_io (s_blank$, r_name, trim(line))

! The exit light, from whichever row holds the run's last state.

pow = 0;  ene = 0;  bun = 0
if (stats%irec > 0) then
  pow = sum(stats%f_power(:, stats%irec));  ene = sum(stats%f_energy(:, stats%irec))
  bun = sum(stats%bunching(:, stats%irec)) / run%nslice
elseif (stats%iend > 0) then
  pow = sum(stats%e_f_power(:, stats%iend));  ene = sum(stats%e_f_energy(:, stats%iend))
  bun = sum(stats%e_bunching(:, stats%iend)) / run%nslice
endif
write (line, '(5a, f6.4)') ' Exit        power ', trim(adjustl(fel_si_str(pow, 'W'))), &
      ', pulse energy ', trim(adjustl(fel_si_str(ene, 'J'))), ', <|b|> ', bun
call out_io (s_blank$, r_name, trim(line))

! The file list is by EXISTENCE, not by replaying which switches were on: a candidate
! name that is there gets listed with its size, and the naming rules stay in the one
! place that owns them (the writers).

n_listed = 0
call note_file (trim(run%global%out_root) // '-final.beam.h5')
do ih = 1, run%n_harm
  hsuf = ''
  if (run%ffield(ih)%harm /= 1) write (hsuf, '(a, i0)') '-h', run%ffield(ih)%harm
  call note_file (trim(run%global%out_root) // '-final' // trim(hsuf) // '.wf.h5')
  call note_file (trim(run%global%out_root) // '-escaped' // trim(hsuf) // '.fld.h5')
  call note_file (trim(run%global%out_root) // '-pulse' // trim(hsuf) // '.fld.h5')
enddo
call note_file (trim(run%global%out_root) // '.stats.h5')
call note_file (trim(run%global%out_root) // '.diag.txt')
call note_file (trim(run%global%out_root) // '.ledger.txt')
call note_file (trim(run%global%out_root) // '.import.txt')
call note_file (trim(run%global%out_root) // '.migration.txt')

call out_io (s_blank$, r_name, &
      '================================================================================')

!------------------------------------------------------------------------------
contains

!+
! Subroutine note_file (name)
!
! Routine to list one output file with its size, if it is there. Existence is the test,
! so the footer never claims a file the run did not write.
!-

subroutine note_file (name)

character(*) name
character(300) row
integer(8) fsize
logical there

!

inquire (file = name, exist = there, size = fsize)
if (.not. there) return
n_listed = n_listed + 1
if (n_listed == 1) then
  write (row, '(2a, t62, a)') ' Wrote       ', trim(name), trim(adjustl(fel_si_str(real(fsize, rp), 'B')))
else
  write (row, '(2a, t62, a)') '             ', trim(name), trim(adjustl(fel_si_str(real(fsize, rp), 'B')))
endif
call out_io (s_blank$, r_name, trim(row))

end subroutine note_file

end subroutine fel_write_footer

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_read_openpmd_into_field (run, fname, ihh, err_flag)
!
! Routine to read an openPMD wavefront into field-set entry ihh (the fundamental
! import path). The photon energy must be the fundamental's: a file carrying a
! harmonic in field_file(1) is refused by name.
!
! Input:
!   run       -- fel_run_struct: Run state.
!   fname     -- character(*): openPMD wavefront file to read.
!   ihh       -- integer: Field-set index to fill.
!
! Output:
!   run       -- fel_run_struct: run%ffield(ihh)%wf holds the imported field.
!   err_flag  -- logical: Set True on a read error or wavelength mismatch. False otherwise.
!-

subroutine fel_read_openpmd_into_field (run, fname, ihh, err_flag)

type (fel_run_struct), target :: run
type (fel_field_struct), pointer :: ffield(:)
type (fel_beam_struct), pointer :: fbeam
character(*) fname
integer ihh
real(rp) e_photon
logical rerr, err_flag
character(*), parameter :: r_name = 'fel_read_openpmd_into_field'

!

err_flag = .true.
ffield => run%ffield
fbeam => run%fbeam

call wavefront_read_openpmd (ffield(ihh)%wf, fname, rerr, e_photon)
if (rerr) return
if (abs(ffield(ihh)%wf%wavelength - fbeam%wavelength) > 1e-6_rp * fbeam%wavelength) then
  call out_io (s_error$, r_name, 'THE openPMD FILE IN FIELD_FILE(1) DOES NOT CARRY THE FUNDAMENTAL:', &
               'ITS PHOTONENERGY WAVELENGTH VS THE BEAM: \2es20.12\ ', &
               r_array = [ffield(ihh)%wf%wavelength, fbeam%wavelength])
  err_flag = .true.;  return
endif
ffield(ihh)%wf%dz = fbeam%slice_spacing
ffield(ihh)%wf%wavelength = fbeam%wavelength   ! One wavelength authority (the 1e-12
                                               ! beam/field consistency check's spirit).

err_flag = .false.

end subroutine fel_read_openpmd_into_field

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_import_harmonic_field (run, fname, err_flag)
!
! Routine to import a harmonic field: match the file's photon energy to the field-set
! entry carrying that harmonic (no match is refused by name), require the fundamental's
! grid and window, and keep the walk's wavelength convention (fundamental / h).
!
! Input:
!   run       -- fel_run_struct: Run state with the fundamental already loaded.
!   fname     -- character(*): openPMD wavefront file to read.
!
! Output:
!   run       -- fel_run_struct: The matching run%ffield entry holds the imported field.
!   err_flag  -- logical: Set True on a read error, no matching harmonic, or a grid
!                  mismatch. False otherwise.
!-

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
character(*), parameter :: r_name = 'fel_import_harmonic_field'

!

err_flag = .true.
ffield => run%ffield
fbeam => run%fbeam
wf => run%ffield(1)%wf
n_harm = run%n_harm

call wavefront_read_openpmd (wtmp, fname, rerr, e_photon)
if (rerr) return
e1 = h_planck * c_light / fbeam%wavelength * e_charge

do ihh = 2, n_harm
  if (abs(e_photon - ffield(ihh)%harm * e1) <= 1e-6_rp * ffield(ihh)%harm * e1) exit
enddo
if (ihh > n_harm) then
  call out_io (s_error$, r_name, 'THE PHOTONENERGY OF ' // trim(fname), &
               '(\es13.5\ J) MATCHES NO FIELD OF THIS RUN''S HARMONICS LIST.', r_array = [e_photon])
  err_flag = .true.;  return
endif

if (size(wtmp%Ex,1) /= size(wf%Ex,1) .or. size(wtmp%Ex,3) /= size(wf%Ex,3) .or. &
    abs(wtmp%dx - wf%dx) > 1e-12_rp * wf%dx) then
  call out_io (s_error$, r_name, 'HARMONIC IMPORT ' // trim(fname), &
               'DOES NOT MATCH THE FUNDAMENTAL''S GRID AND WINDOW (ONE TIME WINDOW, ONE GRID).')
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
!------------------------------------------------------------------------------
!+
! Subroutine fel_drain_bank (run, ihh, err_flag)
!
! Routine to drain the slippage bank: each transmitted slice streams to
! <out_root>-escaped.fld.h5 as it leaves (peak memory a handful of grid planes), with
! its wavefront_params (the full 4x4 at bank time: exactly what the analytic
! free-space propagation of pulse statistics needs) and its transmission z. Genesis
! field-file conventions (dfl units) are used so existing tooling reads the file.
! The root datasets land at finalize.
!
! Input:
!   run       -- fel_run_struct: Run state. run%ffield(ihh)%bank holds the slices to drain.
!   ihh       -- integer: Field-set index whose bank is drained.
!
! Output:
!   run       -- fel_run_struct: Banked params appended to run%bank_z / run%bank_pms.
!                  run%n_banked(ihh) and run%esc_id(ihh) updated.
!   err_flag  -- logical: Set True on a write error. False otherwise.
!-

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
  if (berr) return
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
  if (berr) return
  pms%s = z_now
  run%bank_z(run%n_banked(ihh), ihh) = z_now
  run%bank_pms(:, run%n_banked(ihh), ihh) = [pms%centroid, reshape(pms%sigma, [16]), &
                           pms%energy, pms%power, pms%on_axis_intensity, pms%emit_x, pms%emit_y]
  if (two_pol) then             ! The y component's params add to the banked energy.
    call wavefront_params_of_plane (bnk%plane_y(:,:,k), wfl%dx, wfl%wavelength, &
                                    fbeam%slice_spacing, pms, .true., berr)
    if (berr) return
  endif

  write (gname, '(a, i0.6)') 'slice', run%n_banked(ihh)
  call H5Gcreate_f (run%esc_id(ihh), trim(gname), g_id, h5e)
  if (h5e < 0) return
  work = dfl_scale * reshape(real(bnk%plane(:,:,k), rp), [nx*nx])
  call hdf5_write_dataset_real (g_id, 'field-real', work, berr);  if (berr) return
  work = dfl_scale * reshape(aimag(bnk%plane(:,:,k)), [nx*nx])
  call hdf5_write_dataset_real (g_id, 'field-imag', work, berr);  if (berr) return
  if (two_pol) then
    work = dfl_scale * reshape(real(bnk%plane_y(:,:,k), rp), [nx*nx])
    call hdf5_write_dataset_real (g_id, 'field-real-y', work, berr);  if (berr) return
    work = dfl_scale * reshape(aimag(bnk%plane_y(:,:,k)), [nx*nx])
    call hdf5_write_dataset_real (g_id, 'field-imag-y', work, berr);  if (berr) return
    call hdf5_write_dataset_real (g_id, 'wavefront_params_y', &
            [pms%centroid, reshape(pms%sigma, [16]), pms%energy, pms%power, &
             pms%on_axis_intensity, pms%emit_x, pms%emit_y], berr);  if (berr) return
  endif
  call hdf5_write_dataset_real (g_id, 'z_transmit', [z_now], berr);  if (berr) return
  call hdf5_write_dataset_real (g_id, 'wavefront_params', run%bank_pms(:, run%n_banked(ihh), ihh), berr);  if (berr) return
  call H5Gclose_f (g_id, h5e)
enddo

err_flag = .false.

end subroutine fel_drain_bank

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_escaped_file_name (run, ihh) result (fname)
!
! Routine to construct the escaped-field file name of field ihh. The fundamental
! keeps its pre-harmonic name. A harmonic's file carries -h<h>.
!
! Input:
!   run    -- fel_run_struct: Run state (out_root and the field set).
!   ihh    -- integer: Field-set index.
!
! Output:
!   fname  -- character(420): The escaped-field file name.
!-

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
!------------------------------------------------------------------------------
!+
! Subroutine fel_finalize_diagnostics (run, err_flag)
!
! Routine to finalize the run's file diagnostics: the stats file always, and with
! keep_escaped_field also the escaped file's root datasets and the FULL PULSE at the
! exit plane. For the pulse, each banked slice is read back, free-space propagated over
! z_end - z_transmit (transmitted light is fixed information, and undulator vacuum is
! free space for light with no beam under it), and concatenated above the live
! window. Earliest-transmitted light is furthest ahead, so banked slice k lands at
! pulse index nslice + (run%n_banked - k + 1). The caller has already unrotated the
! live window (the final-dump cshift).
!
! Input:
!   run       -- fel_run_struct: Run state at the end of tracking.
!
! Output:
!   run       -- fel_run_struct: Escaped files closed (run%esc_id zeroed).
!   err_flag  -- logical: Set True on a write error. False otherwise.
!-

subroutine fel_finalize_diagnostics (run, err_flag)

type (fel_run_struct), target :: run
type (fel_field_struct), pointer :: ffield(:)
type (fel_beam_struct), pointer :: fbeam
integer nx, h5e, ihh, n_harm, nslice
logical ferr, err_flag, two_pol, keep_escaped_field
character(400) out_root
character(420) fname_h
character(*), parameter :: r_name = 'fel_finalize_diagnostics'

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
if (ferr) return
call fel_write_meta (run, trim(out_root) // '.stats.h5')

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
    call out_io (s_blank$, r_name, '  ' // trim(fel_escaped_file_name(run, ihh)))
  endif

  ! The full pulse at the exit plane. With two live polarizations, one file per
  ! component (Genesis's format holds one). Harmonic pulses carry -h<h>.

  if (two_pol) then
    call fel_write_pulse_file (run, trim(out_root) // '-pulse-x.fld.h5', .false., ihh, ferr)
    if (ferr) return
    call fel_write_pulse_file (run, trim(out_root) // '-pulse-y.fld.h5', .true., ihh, ferr)
    if (ferr) return
    call out_io (s_blank$, r_name, '  ' // trim(out_root) // '-pulse-{x,y}.fld.h5')
  elseif (ffield(ihh)%harm == 1) then
    call fel_write_pulse_file (run, trim(out_root) // '-pulse.fld.h5', .false., ihh, ferr)
    if (ferr) return
    call out_io (s_blank$, r_name, '  ' // trim(out_root) // '-pulse.fld.h5')
  else
    write (fname_h, '(2a, i0, a)') trim(out_root), '-pulse-h', ffield(ihh)%harm, '.fld.h5'
    call fel_write_pulse_file (run, trim(fname_h), .false., ihh, ferr)
    if (ferr) return
    call out_io (s_blank$, r_name, '  ' // trim(fname_h))
  endif
enddo

err_flag = .false.

end subroutine fel_finalize_diagnostics

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_write_pulse_file (run, fname, use_y, ihh, err_flag)
!
! Routine to write one component of the full exit-plane pulse (Genesis's format holds
! one per file): the live window's slices, then each banked slice read back from the
! escaped file and free-space propagated over z_end - z_transmit. Shared by the -x
! and -y files.
!
! Input:
!   run       -- fel_run_struct: Run state.
!   fname     -- character(*): Pulse file to write.
!   use_y     -- logical: If True write the Ey component, otherwise the Ex component.
!   ihh       -- integer: Field-set index.
!
! Output:
!   err_flag  -- logical: Set True on a read or write error. False otherwise.
!-

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

call hdf5_open_file (trim(fname), 'WRITE', p_id, ferr);  if (ferr) return
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
  call hdf5_write_dataset_real (g_id, 'field-real', work, ferr);  if (ferr) return
  if (use_y) then
    work = dfl_scale * reshape(aimag(wfl%Ey(:,:,is_f)), [nx*nx])
  else
    work = dfl_scale * reshape(aimag(wfl%Ex(:,:,is_f)), [nx*nx])
  endif
  call hdf5_write_dataset_real (g_id, 'field-imag', work, ferr);  if (ferr) return
  call H5Gclose_f (g_id, h5e)
enddo

if (run%n_banked(ihh) > 0) then
  if (.not. allocated(wf1%Ex)) then
    allocate (wf1%Ex(nx, nx, 1))
    wf1%dx = wfl%dx;  wf1%dy = wfl%dy;  wf1%dz = fbeam%slice_spacing
  endif
  wf1%wavelength = wfl%wavelength      ! Per field: banked light drifts at ITS wavelength.

  call hdf5_open_file (trim(fel_escaped_file_name(run, ihh)), 'READ', e_id, ferr);  if (ferr) return
  do k = 1, run%n_banked(ihh)
    write (gname, '(a, i0.6)') 'slice', k
    g_id = hdf5_open_group (e_id, trim(gname), ferr, .true.);  if (ferr) return
    call hdf5_read_dataset_real (g_id, trim(dset_r), re_w, ferr, trim(gname));  if (ferr) return
    call hdf5_read_dataset_real (g_id, trim(dset_i), im_w, ferr, trim(gname));  if (ferr) return
    call H5Gclose_f (g_id, h5e)
    wf1%Ex(:,:,1) = reshape(cmplx(re_w, im_w, wf_rp), [nx, nx]) / dfl_scale

    call wavefront_drift (wf1, z_now - run%bank_z(k, ihh), ferr);  if (ferr) return

    write (gname, '(a, i0.6)') 'slice', nslice + (run%n_banked(ihh) - k + 1)
    call H5Gcreate_f (p_id, trim(gname), g_id, h5e)
    work = dfl_scale * reshape(real(wf1%Ex(:,:,1), rp), [nx*nx])
    call hdf5_write_dataset_real (g_id, 'field-real', work, ferr);  if (ferr) return
    work = dfl_scale * reshape(aimag(wf1%Ex(:,:,1)), [nx*nx])
    call hdf5_write_dataset_real (g_id, 'field-imag', work, ferr);  if (ferr) return
    call H5Gclose_f (g_id, h5e)
  enddo
  call H5Fclose_f (e_id, h5e)
endif

call H5Fclose_f (p_id, h5e)

err_flag = .false.

end subroutine fel_write_pulse_file

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_write_wake_block (run, z)
!
! Routine to write one block of per-slice wake eloss, z-stamped, to the run's open
! wake file. Written at the hoisted update and at every migration-stride recompute.
! The energy-bookkeeping and stale-wake checks parse these blocks.
!
! Input:
!   run  -- fel_run_struct: Run state. run%coll%wake%eloss holds the current losses.
!   z    -- real(rp): The z position stamped on the block [m].
!
! Output:
!   run  -- fel_run_struct: One block appended to the file on unit run%iu_wake.
!-

subroutine fel_write_wake_block (run, z)

type (fel_run_struct), target :: run

real(rp) z
integer is_w

write (run%iu_wake, '(a, es22.14)') '# z = ', z
do is_w = 1, run%nslice
  write (run%iu_wake, '(i8, es24.16)') is_w, run%coll%wake%eloss(is_w)
enddo

end subroutine fel_write_wake_block

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_write_meta (run, stats_file)
!
! Routine to write the provenance group into the stats file (Genesis parity: its Meta
! embeds the entire input and lattice plus timestamp/user/cwd/version). Meta/ carries
! the RESOLVED input echo (every default explicit, straight from the structs), the
! lattice file's text and name, an ISO timestamp, user, cwd and the Bmad version.
! All are ATTRIBUTES, so every dataset-level identity comparison (thread runs,
! re-entrancy passes) is untouched by run-specific provenance.
!
! Note: This is best effort. A failure warns, never fails the run (HDF5 caps compact
! attributes near 64 kB, and a very large lattice would hit it).
!
! Input:
!   run         -- fel_run_struct: Run state.
!   stats_file  -- character(*): The stats file to reopen and annotate.
!
! Output:
!   None. Failures are reported with a warning and otherwise ignored.
!-

subroutine fel_write_meta (run, stats_file)

type (fel_run_struct), target :: run
character(*) stats_file
character(:), allocatable :: txt
character(24) stamp
character(200) user_name, cwd
character(8) date_s
character(10) time_s
integer(hid_t) f_id, g_id
integer h5e
logical merr
character(*), parameter :: r_name = 'fel_write_meta'

!

call hdf5_open_file (stats_file, 'APPEND', f_id, merr)
if (merr) then
  call out_io (s_warn$, r_name, 'Could not reopen the stats file for Meta/.')
  return
endif
call H5Gcreate_f (f_id, 'Meta', g_id, h5e)
if (h5e < 0) then
  call H5Fclose_f (f_id, h5e)
  return
endif

call resolved_input_text (run, txt)
call hdf5_write_attribute_string (g_id, 'input_echo', txt, merr)
if (merr) call out_io (s_warn$, r_name, 'Meta/input_echo did not fit an attribute.')

call file_text (trim(run%lat_file), txt)
call hdf5_write_attribute_string (g_id, 'lattice_text', txt, merr)
if (merr) call out_io (s_warn$, r_name, 'Meta/lattice_text did not fit an attribute.')
call hdf5_write_attribute_string (g_id, 'lattice_file', trim(run%lat_file), merr)

call date_and_time (date_s, time_s)
stamp = date_s(1:4) // '-' // date_s(5:6) // '-' // date_s(7:8) // 'T' // &
        time_s(1:2) // ':' // time_s(3:4) // ':' // time_s(5:6)
call hdf5_write_attribute_string (g_id, 'timestamp', stamp, merr)
call get_environment_variable ('USER', user_name)
call hdf5_write_attribute_string (g_id, 'user', trim(user_name), merr)
call get_environment_variable ('PWD', cwd)
call hdf5_write_attribute_string (g_id, 'cwd', trim(cwd), merr)
call hdf5_write_attribute_int (g_id, 'bmad_inc_version', bmad_inc_version$, merr)

call H5Gclose_f (g_id, h5e)
call H5Fclose_f (f_id, h5e)

end subroutine fel_write_meta

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine resolved_input_text (run, txt)
!
! Routine to render the resolved input echo as one string. Goes through
! fel_write_resolved_input and a scratch file since namelist output needs an
! external unit.
!
! Input:
!   run  -- fel_run_struct: Run state.
!
! Output:
!   txt  -- character(:), allocatable: The resolved input echo.
!-

subroutine resolved_input_text (run, txt)

type (fel_run_struct), target :: run
character(:), allocatable :: txt
integer iu

open (newunit = iu, status = 'scratch', action = 'readwrite')
call fel_write_resolved_input (run, iu)
call unit_text (iu, txt)
close (iu)

end subroutine resolved_input_text

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine file_text (file_name, txt)
!
! Routine to read a whole file into one string. An unreadable file yields the
! string '(unreadable)'.
!
! Input:
!   file_name  -- character(*): File to read.
!
! Output:
!   txt        -- character(:), allocatable: The file's text.
!-

subroutine file_text (file_name, txt)

character(*) file_name
character(:), allocatable :: txt
integer iu, ios

open (newunit = iu, file = file_name, status = 'old', action = 'read', iostat = ios)
if (ios /= 0) then
  txt = '(unreadable)'
  return
endif
call unit_text (iu, txt)
close (iu)

end subroutine file_text

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine unit_text (iu, txt)
!
! Routine to read an open unit from the start into one string.
!
! Input:
!   iu   -- integer: Open unit. Rewound before reading.
!
! Output:
!   txt  -- character(:), allocatable: The unit's text.
!-

subroutine unit_text (iu, txt)

integer iu, ios
character(:), allocatable :: txt
character(1000) line

rewind (iu)
txt = ''
do
  read (iu, '(a)', iostat = ios) line
  if (ios /= 0) exit
  txt = txt // trim(line) // new_line('a')
enddo

end subroutine unit_text

end module fel_io_mod
