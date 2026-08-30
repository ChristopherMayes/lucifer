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
! For bp_com%num_lat_files, the parser's tally of how many files it opened. Read in
! fel_write_meta, which is the only place this module reaches into the parser.
use bmad_parser_struct, only: bp_com
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
! EXT_Wavefront (.wf.h5), both polarizations as components of one mesh record. Genesis
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

call fel_stats_exit_light (stats, pow, ene, bun)
write (line, '(5a, f6.4)') ' Exit        power ', trim(adjustl(fel_si_str(pow, 'W'))), &
      ', pulse energy ', trim(adjustl(fel_si_str(ene, 'J'))), ', <|b|> ', bun
call out_io (s_blank$, r_name, trim(line))

! The file list is by existence, not by replaying which switches were on: a candidate
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
! keep_escaped_field also the escaped file's root datasets and the full pulse at the
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
type (fel_stats_params_struct) sprm
integer nx, h5e, ihh, n_harm, nslice, is
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

! The parameters the file states as data. Assembled here because this is where the
! whole run is visible; fel_stats_mod knows the accumulator and nothing else.

sprm%slice_spacing = fbeam%slice_spacing
sprm%species = 'electron'
sprm%beta0 = fel_p0_mc(fbeam) / sqrt(fel_p0_mc(fbeam)**2 + 1)

call fel_stats_write (run%stats, sprm, trim(out_root) // '.stats.h5', ferr)
call fel_write_params (run, trim(out_root) // '.stats.h5')
if (ferr) return
call fel_write_lattice (run, trim(out_root) // '.stats.h5')
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
  wf1%wavelength = wfl%wavelength      ! Per field: banked light drifts at its wavelength.

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
! Routine to write the provenance group into the stats file: the resolved input echo
! (every default explicit, straight from the structs), the top-level lattice file's name
! and text, how many files the parser opened, an ISO timestamp and the Bmad version.
!
! Datasets, not attributes. HDF5 caps a single attribute at 64 kB and the largest that
! writes here is 65495 bytes, while a scalar string dataset took 3 MB. An echoed namelist
! runs to 12 kB and a lattice text to 37 kB on a real lattice, so the old attributes were
! at half the cap already, and the failure path was a warning: the provenance would have
! gone missing silently on a lattice of no unusual size. The cost is that
! meta/ no longer escapes dataset-level identity comparisons for free, since input_echo
! carries out_root and two runs differ there. The harness excludes meta/ by name instead.
!
! It is not a reproducibility record, and it no longer claims to be. lattice_source is the
! Top-level file only: Bmad's "call, file =" pulls in more, and every wrapper lattice in
! this tree records a call statement while the lattice it calls is absent.
! n_lattice_files says how many files the parser actually opened, so a reader can see at
! a glance whether the text is the whole story. What the file offers for reproduction is
! the lattice/ table (every tracked element with the values the physics used) beside the
! input echo. Serializing the lattice is not an option: write_bmad_lattice_file inlines a
! grid_field as ASCII under one_file$ and writes sibling binary files otherwise, so no
! output_form is both complete and bounded, and Bmad has no HDF5 lattice format to
! borrow. lat%creation_hash is not a substitute either, hashing inode and size, which
! makes it machine specific rather than a fingerprint of content.
!
! A stats file is meant to travel, so nothing here identifies a person by default. The
! timestamp and the Bmad version identify the run. The user name and working directory
! go in only under global%record_environment, and the file records a lattice basenAME.
! What a user types into the namelist is echoed as typed, so an absolute path there is
! still an absolute path in the file. Genesis records user and cwd always. Parity is not
! a reason to leak.
!
! Note: This is best effort. A failure warns, never fails the run.
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
character(400) base, path
character(8) date_s
character(10) time_s
integer(hid_t) f_id, g_id
integer h5e, nfile, ix
logical merr
character(*), parameter :: r_name = 'fel_write_meta'

!

call hdf5_open_file (stats_file, 'APPEND', f_id, merr)
if (merr) then
  call out_io (s_warn$, r_name, 'Could not reopen the stats file for meta/.')
  return
endif
call H5Gcreate_f (f_id, 'meta', g_id, h5e)
if (h5e < 0) then
  call H5Fclose_f (f_id, h5e)
  return
endif
call H5LTset_attribute_string_f (f_id, 'meta', 'kind', 'provenance', h5e)
call H5LTset_attribute_string_f (f_id, 'meta', 'description', &
      'How this file was made. Which lattice, not the lattice itself.', h5e)

merr = .false.
call resolved_input_text (run, txt)
call fel_h5_text (g_id, 'input_echo', 'input echo', &
      'The run''s input with every default made explicit, written from the structs. ' // &
      'File names appear as the user typed them.', txt, merr)

! The lattice: which one, and how much of it is here. bmad_parser opens the top-level
! file plus every file its "call, file =" statements reach, and file_text reads one of
! them, so a wrapper lattice records a call statement and drops what it calls. The count
! is the honesty signal. It comes from the parser's own tally, which outlives the parse
! because parser_end_stuff deallocates only the array of names, and zero means the tally
! was not available rather than that no file was read.

ix = splitfilename(run%lat_file, path, base)
call fel_h5_text (g_id, 'lattice_file', 'lattice file', &
      'Base name of the top-level lattice file. The directory is deliberately absent: ' // &
      'a stats file travels, and a path is machine local.', base, merr)

call file_text (trim(run%lat_file), txt)
call fel_h5_text (g_id, 'lattice_source', 'lattice source', &
      'Text of the top-level lattice file only. A file reached by "call, file =" is ' // &
      'NOT here, so this is not a reproducibility record. See n_lattice_files, and ' // &
      'read lattice/ for the tracked line as the physics used it.', txt, merr)

nfile = max(bp_com%num_lat_files, 0)
call fel_h5_int (g_id, 'n_lattice_files', '1', 'lattice files', &
      'How many files the parser opened, so one means lattice_source is the whole ' // &
      'story. Zero means the parser did not report a count.', '', nfile, merr)

call date_and_time (date_s, time_s)
stamp = date_s(1:4) // '-' // date_s(5:6) // '-' // date_s(7:8) // 'T' // &
        time_s(1:2) // ':' // time_s(3:4) // ':' // time_s(5:6)
call fel_h5_text (g_id, 'timestamp', 'timestamp', &
      'Local time the run finished, ISO 8601 without a zone.', stamp, merr)
call fel_h5_int (g_id, 'bmad_inc_version', '1', 'Bmad version', &
      'The bmad_inc_version$ the run was built against.', '', bmad_inc_version$, merr)

! The environment, only when asked for: these identify a person and a machine, and they
! are of no use to a reader elsewhere.

if (run%global%record_environment) then
  call get_environment_variable ('USER', user_name)
  call fel_h5_text (g_id, 'user', 'user', &
        'Who ran it, from $USER. Present because global%record_environment was set.', &
        user_name, merr)
  call get_environment_variable ('PWD', cwd)
  call fel_h5_text (g_id, 'cwd', 'directory', &
        'Where it ran, from $PWD. Present because global%record_environment was set.', &
        cwd, merr)
endif

if (merr) call out_io (s_warn$, r_name, 'Could not write all of meta/.')

call H5Gclose_f (g_id, h5e)
call H5Fclose_f (f_id, h5e)

end subroutine fel_write_meta

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_write_params (run, stats_file)
!
! Routine to write the input tree into the stats file: params/, one subgroup per input
! structure the program honors, each carrying @struct and holding the resolved value of
! every honored component. What the user controlled, as data, so two runs diff by their
! inputs and no reader needs a defaults table. What the run produced is run/, written by
! fel_stats_write, and the literal text stays in meta/input_echo.
!
! Components the program refuses are not here: recording a knob that did nothing would
! lie about the run (the quiet start's honored beam_init set is doc/input-reference.md's
! table). List-valued inputs (the dump locators, the harmonic set, the field files) are
! string or integer array attributes on their subgroup rather than datasets, since a
! dataset needs an axis and a list of inputs has none worth defining. A blank file-name
! input is omitted: absence means unset. bmad_com and space_charge_com are written
! whole, since the namelist exposes them whole and all of Bmad's tracking reads them.
!
! Best effort, like fel_write_meta: a failure warns and never fails the run.
!
! Input:
!   run         -- fel_run_struct: Run state, after setup.
!   stats_file  -- character(*): The stats file to reopen and annotate.
!-

subroutine fel_write_params (run, stats_file)

type (fel_run_struct), target :: run
character(*) stats_file
integer(hid_t) f_id, p_id, g_id
integer h5e, i, n
logical merr
character(60), allocatable :: locs(:)
character(400), allocatable :: files(:)
character(*), parameter :: r_name = 'fel_write_params'

!

call hdf5_open_file (stats_file, 'APPEND', f_id, merr)
if (merr) then
  call out_io (s_warn$, r_name, 'Could not reopen the stats file for params/.')
  return
endif
call H5Gcreate_f (f_id, 'params', p_id, h5e)
if (h5e < 0) then
  call H5Fclose_f (f_id, h5e)
  return
endif
merr = .false.
call H5LTset_attribute_string_f (f_id, 'params', 'kind', 'input', h5e)
call H5LTset_attribute_string_f (f_id, 'params', 'description', &
      'What the user controlled: one subgroup per honored input structure, every ' // &
      'honored component resolved after defaults.', h5e)

! ------------------------------------------------------------ global.

! out_root is NOT here, deliberately: it is the run's own name rather than a physics
! input. Every output file already carries it, the echo records it as typed, and as
! data it would make two otherwise identical runs compare different, which is exactly
! what the thread-identity check must not see.

call sub_open ('global', 'fel_global_struct', 'The run switches (&fel_params global%).')
call fel_h5_str (g_id, 'interlude_model', 'interludes', &
      'Interlude transport: bmad (the seam) or genesis (transcribed).', '', &
      [run%global%interlude_model], merr)
call fel_h5_str (g_id, 'source_model', 'source', &
      'The FEL source model: deposit or coherent.', '', [run%global%source_model], merr)
call fel_h5_text (g_id, 'track_start', 'track start', &
      'Element locator bounding the walk. Blank means the whole line.', &
      run%global%track_start, merr)
call fel_h5_text (g_id, 'track_end', 'track end', &
      'Element locator bounding the walk. Blank means the whole line.', &
      run%global%track_end, merr)
call fel_h5_real (g_id, 'comb_ds_save', 'm', 'comb', &
      'Minimum z advance between per-record rows. 0 every record, negative none.', '', &
      run%global%comb_ds_save, merr)
call fel_h5_int (g_id, 'ran_seed', '1', 'seed', &
      'The one RNG seed: generation, import and noise.', '', run%global%ran_seed, merr)
call fel_h5_flag (g_id, 'write_diag', 'write diag', &
      'Write the Genesis-comparison text diag file.', run%global%write_diag, merr)
call fel_h5_flag (g_id, 'write_initial', 'write initial', &
      'Dump the initial state before tracking.', run%global%write_initial, merr)
call fel_h5_flag (g_id, 'load_only', 'load only', &
      'Build the initial state, dump it, stop.', run%global%load_only, merr)
call fel_h5_flag (g_id, 'keep_escaped_field', 'keep escaped', &
      'Bank the field slices slippage transmits out of the window.', &
      run%global%keep_escaped_field, merr)
call fel_h5_flag (g_id, 'migrate', 'migrate', &
      'Move particles between slices when their phase leaves the slice window.', &
      run%global%migrate, merr)
call fel_h5_flag (g_id, 'migrate_check', 'migrate check', &
      'Verify phase continuity at each migration.', run%global%migrate_check, merr)
call fel_h5_flag (g_id, 'reference_run', 'reference run', &
      'No FEL interaction: Bmad tracks everything.', run%global%reference_run, merr)
call fel_h5_flag (g_id, 'record_environment', 'record environment', &
      'Record the user name and working directory in meta/.', &
      run%global%record_environment, merr)
call loc_note ('dump_beam_at', run%global%dump_beam_at)
call loc_note ('dump_field_at', run%global%dump_field_at)
call H5Gclose_f (g_id, h5e)

! ------------------------------------------------------------ beam_init, honored set.

call sub_open ('beam_init', 'beam_init_struct', &
      'The bunch description (&fel_beam_init beam_init%), the quiet start''s honored set.')
call fel_h5_int (g_id, 'n_particle', '1', 'particles per slice', &
      'Macroparticles per slice (the import path reads it as bunch particles).', '', &
      run%beam_init%n_particle, merr)
call fel_h5_real (g_id, 'bunch_charge', 'C', 'charge', &
      'Total bunch charge. The current is derived from it, never input.', '', &
      run%beam_init%bunch_charge, merr)
call fel_h5_real (g_id, 'a_norm_emit', 'm rad', 'a norm emit', &
      'Normalized a-mode emittance.', '', run%beam_init%a_norm_emit, merr)
call fel_h5_real (g_id, 'b_norm_emit', 'm rad', 'b norm emit', &
      'Normalized b-mode emittance.', '', run%beam_init%b_norm_emit, merr)
call fel_h5_real (g_id, 'sig_pz', '1', 'sig_pz', &
      'Fractional momentum spread dP/P0.', '', run%beam_init%sig_pz, merr)
call fel_h5_real (g_id, 'sig_z', 'm', 'sig_z', &
      'Bunch length. Zero is the steady state for a one-slice window.', '', &
      run%beam_init%sig_z, merr)
call fel_h5_str (g_id, 'distribution_type', 'distribution', &
      'Bmad''s per-plane distribution names, the z entry selecting the current ' // &
      'profile: RAN_GAUSS is Gaussian, GRID is flat.', 'plane', &
      run%beam_init%distribution_type, merr)
call fel_h5_real (g_id, 'grid_z_min', 'm', 'grid z min', &
      'Flat-profile z extent, beam_init%grid(3)%x_min. Read under GRID only.', '', &
      run%beam_init%grid(3)%x_min, merr)
call fel_h5_real (g_id, 'grid_z_max', 'm', 'grid z max', &
      'Flat-profile z extent, beam_init%grid(3)%x_max. Read under GRID only.', '', &
      run%beam_init%grid(3)%x_max, merr)
call H5Gclose_f (g_id, h5e)

! ------------------------------------------------------------ the beam-side scalars.

call sub_open ('beam_param', 'fel_beam_init_param_struct', &
      'The beam-side scalars of &fel_beam_init beside beam_init and imp.')
call fel_h5_int (g_id, 'nbins', '1', 'beamlet size', &
      'Beamlet size of the quiet start (quiet below nbins).', '', run%bparam%nbins, merr)
call fel_h5_flag (g_id, 'shotnoise', 'shot noise', &
      'Impose physical (Fawley) shot noise on the phases.', run%bparam%shotnoise, merr)
call fel_h5_flag (g_id, 'use_beam_init', 'use beam_init', &
      'Generate the bunch from beam_init, then import it.', run%bparam%use_beam_init, merr)
call fel_h5_flag (g_id, 'split_weights', 'split weights', &
      'Check instrument: coincident w/3 + 2w/3 copies after loading.', &
      run%bparam%split_weights, merr)
call fel_h5_flag (g_id, 'swap_beam_xy', 'swap xy', &
      'Check instrument: swap (x,px) and (y,py) after generation.', &
      run%bparam%swap_beam_xy, merr)
call fel_h5_flag (g_id, 'gen_test_weights', 'test weights', &
      'Check instrument: alternate beamlet weights 0.25x/1.75x.', &
      run%bparam%gen_test_weights, merr)
call fel_h5_flag (g_id, 'imp_split_weights', 'import split', &
      'Check instrument: split-weight copies before the import resample.', &
      run%bparam%imp_split_weights, merr)
call file_note ('beam_file', run%bparam%beam_file, 'Genesis .par.h5 dump to import.')
call file_note ('dist_file', run%bparam%dist_file, 'openPMD-beamphysics particle file to import.')
call file_note ('write_dist_file', run%bparam%write_dist_file, 'Write the bunch as a Genesis DISTRIBUTION file.')
call file_note ('write_opmd_file', run%bparam%write_opmd_file, 'Write the bunch as openPMD-beamphysics.')
call H5Gclose_f (g_id, h5e)

! ------------------------------------------------------------ the import resampler.

call sub_open ('imp', 'fel_import_param_struct', &
      'The importdistribution resampler knobs (&fel_beam_init imp%).')
call fel_h5_real (g_id, 'slicewidth', '1', 'slice width', &
      'Sampling window over bunch length (Genesis''s slicewidth).', '', run%imp%slicewidth, merr)
call fel_h5_int (g_id, 'npart', '1', 'particles', &
      'Macroparticles per slice after resampling.', '', run%imp%npart, merr)
call fel_h5_int (g_id, 'nbins', '1', 'beamlet size', &
      'Beamlet size of the resample.', '', run%imp%nbins, merr)
call fel_h5_int (g_id, 'nslice', '1', 'slices', &
      'Slice count. Zero means round(bunch length / spacing), Genesis''s rule.', '', &
      run%imp%nslice, merr)
call H5Gclose_f (g_id, h5e)

! ------------------------------------------------------------ the radiation start.

call sub_open ('wavefront_init', 'wavefront_init_struct', &
      'The radiation starting condition (&fel_wavefront_init wavefront_init%).')
call fel_h5_real (g_id, 'lambda0', 'm', 'lambda0', &
      'Fundamental radiation wavelength.', '', run%winit%lambda0, merr)
call fel_h5_real (g_id, 'window_length', 'm', 'window length', &
      'Time window. Zero means derived from the bunch.', '', run%winit%window_length, merr)
call fel_h5_int (g_id, 'window_sample', '1', 'sample', &
      'Slice spacing in wavelengths (Genesis''s sample), and so the number of ' // &
      'undulator periods of slippage per slice. An integer, which is what lets the ' // &
      'field record rotate by one index with no interpolation.', '', &
      run%winit%window_sample, merr)
call fel_h5_int (g_id, 'grid_n_pts', '1', 'grid points', &
      'Transverse grid points per side.', '', run%winit%grid_n_pts, merr)
call fel_h5_real (g_id, 'grid_half_width', 'm', 'grid half width', &
      'Transverse grid half width.', '', run%winit%grid_half_width, merr)
call fel_h5_real (g_id, 'seed_power', 'W', 'seed power', &
      'Gaussian seed power. Zero is a dark start.', '', run%winit%seed_power, merr)
call fel_h5_real (g_id, 'seed_waist_size', 'm', 'seed waist', &
      'Seed intensity 1/e^2 radius.', '', run%winit%seed_waist_size, merr)
call fel_h5_str (g_id, 'seed_polarization', 'polarization', &
      'Seed polarization, x or y.', '', [run%winit%seed_polarization], merr)
n = count(run%winit%harmonics > 0)
call H5LTset_attribute_int_f (g_id, '.', 'harmonics', run%winit%harmonics(1:n), &
                              int(n, size_t), h5e)
merr = merr .or. (h5e < 0)
allocate (files(count(run%field_file /= '')))
n = 0
do i = 1, size(run%field_file)
  if (run%field_file(i) /= '') then
    n = n + 1;  files(n) = run%field_file(i)
  endif
enddo
if (n > 0) then
  call hdf5_write_attribute_string_rank1 (g_id, 'field_file', files(1:n), merr)
endif
deallocate (files)
call H5Gclose_f (g_id, h5e)

! ------------------------------------------------------------ the chamber wake.

call sub_open ('wake', 'fel_wake_init_struct', &
      'The chamber-wake description (&fel_params wake%), Genesis &wake names.')
call fel_h5_flag (g_id, 'on', 'wake on', 'Any chamber wake at all.', run%wake_init%on, merr)
call fel_h5_real (g_id, 'loss', 'eV/m', 'loss', 'External loss.', '', run%wake_init%loss, merr)
call fel_h5_real (g_id, 'radius', 'm', 'radius', &
      'Chamber radius, or half gap if flat.', '', run%wake_init%radius, merr)
call fel_h5_real (g_id, 'conductivity', '1/(Ohm m)', 'conductivity', &
      'DC conductivity. Zero means no resistive wake.', '', run%wake_init%conductivity, merr)
call fel_h5_real (g_id, 'relaxation', 'm', 'relaxation', &
      'AC relaxation distance c*tau.', '', run%wake_init%relaxation, merr)
call fel_h5_flag (g_id, 'roundpipe', 'round pipe', &
      'Round chamber. False is flat, parallel plates.', run%wake_init%roundpipe, merr)
call fel_h5_str (g_id, 'material', 'material', &
      'CU or AL shortcut for conductivity and relaxation. Blank when set directly.', &
      '', [run%wake_init%material], merr)
call fel_h5_real (g_id, 'gap', 'm', 'gap', &
      'Undulator gap. Zero means no geometric wake.', '', run%wake_init%gap, merr)
call fel_h5_real (g_id, 'lgap', 'm', 'gap period', 'Period of the gaps.', '', &
      run%wake_init%lgap, merr)
call fel_h5_real (g_id, 'hrough', 'm', 'roughness', &
      'Roughness amplitude. Zero means no roughness wake.', '', run%wake_init%hrough, merr)
call fel_h5_real (g_id, 'lrough', 'm', 'roughness period', 'Roughness period.', '', &
      run%wake_init%lrough, merr)
call file_note ('write_kernels', run%wake_init%write_kernels, &
      'Check instrument: export the transcribed kernels.')
call H5Gclose_f (g_id, h5e)

! ------------------------------------------------------------ space charge.

call sub_open ('sc', 'fel_efield_struct', &
      'The FEL space-charge description (&fel_params sc%).')
call fel_h5_flag (g_id, 'on', 'space charge on', 'Any space charge at all.', &
      run%sc_init%on, merr)
call fel_h5_real (g_id, 'rmax', 'm', 'rmax', &
      'Radial grid extent scale. Grows adaptively as Genesis''s.', '', run%sc_init%rmax, merr)
call fel_h5_int (g_id, 'ngrid', '1', 'radial points', 'Radial grid points.', '', &
      run%sc_init%ngrid, merr)
call fel_h5_int (g_id, 'nz', '1', 'harmonics', &
      'Longitudinal harmonics. Zero disables the short-range solve.', '', run%sc_init%nz, merr)
call fel_h5_int (g_id, 'nphi', '1', 'azimuthal modes', &
      'Azimuthal modes m = -nphi..nphi.', '', run%sc_init%nphi, merr)
call fel_h5_flag (g_id, 'longrange', 'long range', &
      'The whole-window longESC term.', run%sc_init%longrange, merr)
call H5Gclose_f (g_id, h5e)

! ------------------------------------------------------------ bmad_com, whole.

call write_bmad_com ()
call write_space_charge_com ()

if (merr) call out_io (s_warn$, r_name, 'Could not write all of params/.')
call H5Gclose_f (p_id, h5e)
call H5Fclose_f (f_id, h5e)

!------------------------------------------------------------------------------
contains

!+
! Subroutine sub_open (name, struct_name, descrip)
!
! Routine to open one input subgroup with its three attributes: kind, @struct naming
! the Fortran type verbatim, and the description.
!-

subroutine sub_open (name, struct_name, descrip)

integer h5e2
character(*) name, struct_name, descrip

!

call H5Gcreate_f (p_id, name, g_id, h5e2)
merr = merr .or. (h5e2 < 0)
call H5LTset_attribute_string_f (p_id, name, 'kind', 'input', h5e2)
call H5LTset_attribute_string_f (p_id, name, 'struct', struct_name, h5e2)
call H5LTset_attribute_string_f (p_id, name, 'description', descrip, h5e2)
merr = merr .or. (h5e2 < 0)

end subroutine sub_open

!+
! Subroutine loc_note (attrib, list)
!
! Routine to record a locator list as a string-array attribute of the nonblank
! entries. An attribute rather than a dataset, since a dataset needs an axis and a
! list of inputs has none worth defining. Omitted when empty: absence means unset.
!-

subroutine loc_note (attrib, list)

integer nn, ii
character(*) attrib, list(:)

!

nn = 0
if (allocated(locs)) deallocate (locs)
allocate (locs(count(list /= '')))
do ii = 1, size(list)
  if (list(ii) /= '') then
    nn = nn + 1;  locs(nn) = list(ii)
  endif
enddo
if (nn > 0) call hdf5_write_attribute_string_rank1 (g_id, attrib, locs(1:nn), merr)

end subroutine loc_note

!+
! Subroutine file_note (name, val, descrip)
!
! Routine to record one file-name input as a text dataset, as the user typed it.
! A blank one is omitted: absence means unset.
!-

subroutine file_note (name, val, descrip)

character(*) name, val, descrip

!

if (val == '') return
call fel_h5_text (g_id, name, name, descrip // ' As the user typed it.', val, merr)

end subroutine file_note

!+
! Subroutine write_bmad_com ()
!
! Routine to write Bmad's bmad_com whole: the namelist exposes every component and all
! of Bmad's tracking reads them, so every one is honored. Descriptions follow the
! struct's own comments.
!-

subroutine write_bmad_com ()

!

call sub_open ('bmad_com', 'bmad_common_struct', &
      'Bmad''s global tracking switches (&fel_params bmad_com%), written whole.')
call fel_h5_real (g_id, 'max_aperture_limit', 'm', 'max aperture', 'Max aperture.', '', &
      bmad_com%max_aperture_limit, merr)
call fel_h5_real (g_id, 'd_orb', 'm,1,m,1,m,1', 'orbit deltas', &
      'Orbit deltas for the mat6-via-tracking calc, per phase-space coordinate.', &
      'bmad', bmad_com%d_orb, merr)
call H5LTset_attribute_string_f (g_id, 'd_orb', 'unit_of_axis', 'bmad', h5e)
merr = merr .or. (h5e < 0)
call fel_h5_dset_attr_int (g_id, 'd_orb', 'unit_power', 1, merr)
call fel_h5_real (g_id, 'default_ds_step', 'm', 'default step', &
      'Default integration step for elements without an explicit one.', '', &
      bmad_com%default_ds_step, merr)
call fel_h5_real (g_id, 'significant_length', 'm', 'significant length', &
      'Below this a length is zero.', '', bmad_com%significant_length, merr)
call fel_h5_real (g_id, 'rel_tol_tracking', '1', 'rel tol', &
      'Closed-orbit relative tolerance.', '', bmad_com%rel_tol_tracking, merr)
call fel_h5_real (g_id, 'abs_tol_tracking', '1', 'abs tol', &
      'Closed-orbit absolute tolerance.', '', bmad_com%abs_tol_tracking, merr)
call fel_h5_real (g_id, 'rel_tol_adaptive_tracking', '1', 'RK rel tol', &
      'Runge-Kutta tracking relative tolerance.', '', bmad_com%rel_tol_adaptive_tracking, merr)
call fel_h5_real (g_id, 'abs_tol_adaptive_tracking', '1', 'RK abs tol', &
      'Runge-Kutta tracking absolute tolerance.', '', bmad_com%abs_tol_adaptive_tracking, merr)
call fel_h5_real (g_id, 'init_ds_adaptive_tracking', 'm', 'RK initial step', &
      'Initial adaptive step size.', '', bmad_com%init_ds_adaptive_tracking, merr)
call fel_h5_real (g_id, 'min_ds_adaptive_tracking', 'm', 'RK min step', &
      'Minimum adaptive step size.', '', bmad_com%min_ds_adaptive_tracking, merr)
call fel_h5_real (g_id, 'fatal_ds_adaptive_tracking', 'm', 'RK fatal step', &
      'Below this step the particle is lost.', '', bmad_com%fatal_ds_adaptive_tracking, merr)
call fel_h5_real (g_id, 'autoscale_amp_abs_tol', 'eV', 'autoscale abs', &
      'Autoscale absolute amplitude tolerance.', '', bmad_com%autoscale_amp_abs_tol, merr)
call fel_h5_real (g_id, 'autoscale_amp_rel_tol', '1', 'autoscale rel', &
      'Autoscale relative amplitude tolerance.', '', bmad_com%autoscale_amp_rel_tol, merr)
call fel_h5_real (g_id, 'autoscale_phase_tol', 'rad', 'autoscale phase', &
      'Autoscale phase tolerance.', '', bmad_com%autoscale_phase_tol, merr)
call fel_h5_real (g_id, 'electric_dipole_moment', '1', 'EDM', &
      'Particle electric dipole moment.', '', bmad_com%electric_dipole_moment, merr)
call fel_h5_real (g_id, 'synch_rad_scale', '1', 'synch rad scale', &
      'Synchrotron radiation kick scale. One is physical, zero is off.', '', &
      bmad_com%synch_rad_scale, merr)
call fel_h5_real (g_id, 'sad_eps_scale', '1', 'sad eps', &
      'sad_mult step length scale.', '', bmad_com%sad_eps_scale, merr)
call fel_h5_real (g_id, 'sad_amp_max', '1', 'sad amp max', &
      'sad_mult step length amplitude cap.', '', bmad_com%sad_amp_max, merr)
call fel_h5_int (g_id, 'sad_n_div_max', '1', 'sad divisions', &
      'sad_mult maximum divisions.', '', bmad_com%sad_n_div_max, merr)
call fel_h5_int (g_id, 'taylor_order', '1', 'taylor order', &
      'Taylor order. Zero means PTC''s saved default.', '', bmad_com%taylor_order, merr)
call fel_h5_int (g_id, 'runge_kutta_order', '1', 'RK order', &
      'Runge-Kutta order.', '', bmad_com%runge_kutta_order, merr)
call fel_h5_int (g_id, 'default_integ_order', '1', 'PTC order', &
      'PTC integration order.', '', bmad_com%default_integ_order, merr)
call fel_h5_int (g_id, 'max_num_runge_kutta_step', '1', 'RK step cap', &
      'Maximum RK steps before the particle counts as lost.', '', &
      bmad_com%max_num_runge_kutta_step, merr)
call fel_h5_flag (g_id, 'rf_phase_below_transition_ref', 'below transition', &
      'Autoscale uses the below-transition stable point.', &
      bmad_com%rf_phase_below_transition_ref, merr)
call fel_h5_flag (g_id, 'sr_wakes_on', 'sr wakes', 'Short-range wakefields.', &
      bmad_com%sr_wakes_on, merr)
call fel_h5_flag (g_id, 'lr_wakes_on', 'lr wakes', 'Long-range wakefields.', &
      bmad_com%lr_wakes_on, merr)
call fel_h5_flag (g_id, 'auto_bookkeeper', 'auto bookkeeper', &
      'Deprecated and no longer used.', bmad_com%auto_bookkeeper, merr)
call fel_h5_flag (g_id, 'high_energy_space_charge_on', 'HE space charge', &
      'High-energy space charge.', bmad_com%high_energy_space_charge_on, merr)
call fel_h5_flag (g_id, 'high_energy_space_charge_linear', 'HE SC linear', &
      'High-energy space charge, linear form.', &
      bmad_com%high_energy_space_charge_linear, merr)
call fel_h5_flag (g_id, 'csr_and_space_charge_on', 'CSR + SC', &
      'CSR and space charge.', bmad_com%csr_and_space_charge_on, merr)
call fel_h5_flag (g_id, 'spin_tracking_on', 'spin tracking', 'Spin tracking.', &
      bmad_com%spin_tracking_on, merr)
call fel_h5_flag (g_id, 'spin_sokolov_ternov_flipping_on', 'ST flipping', &
      'Sokolov-Ternov spin flipping on emission.', &
      bmad_com%spin_sokolov_ternov_flipping_on, merr)
call fel_h5_flag (g_id, 'radiation_damping_on', 'radiation damping', &
      'Radiation damping.', bmad_com%radiation_damping_on, merr)
call fel_h5_flag (g_id, 'radiation_zero_average', 'damping zero average', &
      'Shift damping to zero on the zero orbit.', bmad_com%radiation_zero_average, merr)
call fel_h5_flag (g_id, 'radiation_fluctuations_on', 'radiation fluctuations', &
      'Radiation fluctuations.', bmad_com%radiation_fluctuations_on, merr)
call fel_h5_flag (g_id, 'conserve_taylor_maps', 'conserve maps', &
      'Bookkeeper may set taylor_map_includes_offsets false.', &
      bmad_com%conserve_taylor_maps, merr)
call fel_h5_flag (g_id, 'absolute_time_tracking', 'absolute time', &
      'Absolute rather than relative time tracking.', &
      bmad_com%absolute_time_tracking, merr)
call fel_h5_flag (g_id, 'absolute_time_ref_shift', 'abs time ref shift', &
      'Apply the reference time shift under absolute time.', &
      bmad_com%absolute_time_ref_shift, merr)
call fel_h5_flag (g_id, 'convert_to_kinetic_momentum', 'kinetic momentum', &
      'Cancel vector-potential kicks in symplectic tracking.', &
      bmad_com%convert_to_kinetic_momentum, merr)
call fel_h5_flag (g_id, 'normalize_twiss', 'normalize twiss', &
      'Normalize the matrix for off-energy twiss.', bmad_com%normalize_twiss, merr)
call fel_h5_flag (g_id, 'aperture_limit_on', 'apertures', &
      'Use apertures in tracking.', bmad_com%aperture_limit_on, merr)
call fel_h5_flag (g_id, 'spin_n0_direction_user_set', 'n0 user set', &
      'User sets n0 for closed branches.', bmad_com%spin_n0_direction_user_set, merr)
call fel_h5_flag (g_id, 'debug', 'debug', 'Code debugging.', bmad_com%debug, merr)
call H5Gclose_f (g_id, h5e)

end subroutine write_bmad_com

!+
! Subroutine write_space_charge_com ()
!
! Routine to write Bmad's space_charge_com whole, same reasoning as bmad_com.
!-

subroutine write_space_charge_com ()

!

call sub_open ('space_charge_com', 'space_charge_common_struct', &
      'Bmad''s space-charge and CSR switches (&fel_params space_charge_com%), whole.')
call fel_h5_real (g_id, 'ds_track_step', 'm', 'CSR step', 'CSR tracking step size.', '', &
      space_charge_com%ds_track_step, merr)
call fel_h5_real (g_id, 'dt_track_step', 's', 'time RK step', &
      'Time Runge-Kutta initial step.', '', space_charge_com%dt_track_step, merr)
call fel_h5_real (g_id, 'cathode_strength_cutoff', '1', 'cathode cutoff', &
      'Cutoff for the cathode field calc.', '', space_charge_com%cathode_strength_cutoff, merr)
call fel_h5_real (g_id, 'rel_tol_tracking', '1', 'rel tol', &
      'Tracking relative tolerance.', '', space_charge_com%rel_tol_tracking, merr)
call fel_h5_real (g_id, 'abs_tol_tracking', '1', 'abs tol', &
      'Tracking absolute tolerance.', '', space_charge_com%abs_tol_tracking, merr)
call fel_h5_real (g_id, 'beam_chamber_height', 'm', 'chamber height', &
      'Used in the shielding calculation.', '', space_charge_com%beam_chamber_height, merr)
call fel_h5_real (g_id, 'lsc_sigma_cutoff', '1', 'LSC cutoff', &
      'Cutoff for the 1D longitudinal SC calc.', '', space_charge_com%lsc_sigma_cutoff, merr)
call fel_h5_real (g_id, 'particle_sigma_cutoff', '1', 'particle cutoff', &
      '3D SC cutoff for far-out particles. Nonpositive means ignore.', '', &
      space_charge_com%particle_sigma_cutoff, merr)
call fel_h5_real (g_id, 'mesh_growth_factor', '1', 'mesh growth', &
      'Fractional padding when growing the SC mesh.', '', &
      space_charge_com%mesh_growth_factor, merr)
call fel_h5_real (g_id, 'mesh_shrink_factor', '1', 'mesh shrink', &
      'Fractional threshold for shrinking the SC mesh.', '', &
      space_charge_com%mesh_shrink_factor, merr)
call fel_h5_int (g_id, 'space_charge_mesh_size', '1', 'SC mesh', &
      'Grid size of the fft_3d space-charge calc, per spatial dimension.', 'plane', &
      space_charge_com%space_charge_mesh_size, merr)
call fel_h5_int (g_id, 'csr3d_mesh_size', '1', 'CSR mesh', &
      'Grid size of the 3D CSR calc, per spatial dimension.', 'plane', &
      space_charge_com%csr3d_mesh_size, merr)
call fel_h5_int (g_id, 'n_bin', '1', 'bins', 'Number of bins used.', '', &
      space_charge_com%n_bin, merr)
call fel_h5_int (g_id, 'particle_bin_span', '1', 'bin span', &
      'Longitudinal particle length over the bin size.', '', &
      space_charge_com%particle_bin_span, merr)
call fel_h5_int (g_id, 'n_shield_images', '1', 'shield images', &
      'Chamber-wall shielding images. Zero means none.', '', &
      space_charge_com%n_shield_images, merr)
call fel_h5_int (g_id, 'sc_min_in_bin', '1', 'min in bin', &
      'Minimum particles per bin for valid sigmas.', '', space_charge_com%sc_min_in_bin, merr)
call fel_h5_flag (g_id, 'lsc_kick_transverse_dependence', 'LSC transverse', &
      'Transverse dependence in the LSC kick.', &
      space_charge_com%lsc_kick_transverse_dependence, merr)
call fel_h5_flag (g_id, 'debug', 'debug', 'Code debugging.', space_charge_com%debug, merr)
call file_note ('diagnostic_output_file', space_charge_com%diagnostic_output_file, &
      'Diagnostic (wake) output file.')
call H5Gclose_f (g_id, h5e)

end subroutine write_space_charge_com

end subroutine fel_write_params

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_write_lattice (run, stats_file)
!
! Routine to write the lattice table into the stats file: one row per tracked element,
! which is what a layout plot needs and what the file has never carried. Genesis writes
! its Lattice/ group as per-step arrays; this is a table instead, joined to the records
! through coords/ix_ele, so nothing is duplicated per record and the two cannot drift.
!
! NOT a lattice serialization. It carries what a plot and a join need, and anything
! that needs the real lattice reads meta/lattice_file or meta/lattice_text, which stay
! the reproducibility record. Bmad has no portable HDF5 lattice format, and a
! statistics file is the wrong place to invent one.
!
! The table's axis is named ele and gets its own coordinate here, coords/ele = 0..n_ele.
! coords/ix_ele is then a variable on the record axis whose values index that one, which
! is the join. Both were called ix_ele before, and one name for two axes of different
! length is a collision every reader has to resolve by hand.
!
! aw is the value the physics used, derived through b_max and l_period with the
! helical or planar factor, not the raw attribute: that is the number a Genesis user
! compares against Lattice/aw, and deriving it here would repeat fel_setup_lattice.
!
! Best effort, like fel_write_meta: a failure warns and never fails the run.
!
! Input:
!   run         -- fel_run_struct: Run state, after setup.
!   stats_file  -- character(*): The stats file to reopen and annotate.
!-

subroutine fel_write_lattice (run, stats_file)

type (fel_run_struct), target :: run
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
character(*) stats_file
character(fel_h5_str_len$), allocatable :: names(:), keys(:)   ! 0:n_ele_track
real(rp), allocatable :: s1(:), s2(:), len_(:), dstep(:), b_max(:), aw(:), l_per(:)
real(rp), allocatable :: ku(:), k1(:), tilt(:), z_off(:)
integer, allocatable :: is_fel(:), helical(:), mode(:)
integer(hid_t) f_id, g_id
integer h5e, ie, ne
logical merr
character(*), parameter :: r_name = 'fel_write_lattice'

!

! The table is indexed BY ix_ele, element 0 (Bmad's beginning element) included, so
! coords/ix_ele indexes it with no offset arithmetic anywhere: lattice/name(ix_ele) is
! the record's element. The first record of a run sits at element 0.

branch => run%lat%branch(0)
ne = branch%n_ele_track

allocate (names(0:ne), keys(0:ne), s1(0:ne), s2(0:ne), len_(0:ne), dstep(0:ne))
allocate (b_max(0:ne), aw(0:ne), l_per(0:ne), ku(0:ne), k1(0:ne), tilt(0:ne), z_off(0:ne))
allocate (is_fel(0:ne), helical(0:ne), mode(0:ne))

do ie = 0, ne
  ele => branch%ele(ie)
  names(ie) = ele%name
  keys(ie) = key_name(ele%key)
  s1(ie) = ele%s_start
  s2(ie) = ele%s
  len_(ie) = ele%value(l$)
  dstep(ie) = ele%value(ds_step$)
  b_max(ie) = ele%value(b_max$)
  l_per(ie) = ele%value(l_period$)
  k1(ie) = ele%value(k1$)
  tilt(ie) = ele%value(tilt_tot$)
  z_off(ie) = ele%value(z_offset_tot$)
  is_fel(ie) = 0
  if (ie > 0) is_fel(ie) = merge(1, 0, run%is_fel(ie))
  aw(ie) = 0;  ku(ie) = 0;  helical(ie) = 0;  mode(ie) = 0
  if (is_fel(ie) == 1) then
    aw(ie) = run%und_of(ie)%aw
    ku(ie) = run%und_of(ie)%ku
    helical(ie) = merge(1, 0, run%und_of(ie)%helical)
    mode(ie) = run%fel_mode(ie)
  endif
enddo

call hdf5_open_file (stats_file, 'APPEND', f_id, merr)
if (merr) then
  call out_io (s_warn$, r_name, 'Could not reopen the stats file for lattice/.')
  return
endif
call H5Gcreate_f (f_id, 'lattice', g_id, h5e)
if (h5e < 0) then
  call H5Fclose_f (f_id, h5e)
  return
endif

merr = .false.
call H5LTset_attribute_string_f (f_id, 'lattice', 'kind', 'table', h5e)
call H5LTset_attribute_string_f (f_id, 'lattice', 'description', &
      'One row per tracked element, indexed by the ele axis.', h5e)

call fel_h5_str (g_id, 'name', 'name', &
      'Element name. Indexed BY the ele axis, element 0 included.', 'ele', names, merr)
call fel_h5_str (g_id, 'key', 'class', 'Bmad element class, as key_name gives it.', &
      'ele', keys, merr)
call fel_h5_real (g_id, 's_start', 'm', 's start', 'Upstream end of the element.', &
      'ele', s1, merr)
call fel_h5_real (g_id, 's_end', 'm', 's end', 'Downstream end of the element.', 'ele', s2, merr)
call fel_h5_real (g_id, 'l', 'm', 'length', 'Element length.', 'ele', len_, merr)
call fel_h5_real (g_id, 'ds_step', 'm', 'ds_step', &
      'Integration step the walk used, which is what sets the record density.', &
      'ele', dstep, merr)
call fel_h5_flag (g_id, 'is_fel', 'is FEL', &
      'One where the element is an FEL segment the FEL step tracked.', 'ele', is_fel, merr)
call fel_h5_int (g_id, 'fel_tracking', '1', 'tracking mode', &
      'Tracking mode of an FEL segment: -1 transcribed Genesis maps, 0 averaged ' // &
      '(the default, Bmad''s own kernel), 1 unaveraged. Zero off an FEL segment.', &
      'ele', mode, merr)
call fel_h5_real (g_id, 'b_max', 'T', 'b_max', 'Peak undulator field, zero elsewhere.', &
      'ele', b_max, merr)
call fel_h5_real (g_id, 'aw', '1', 'aw', &
      'Rms undulator parameter as the physics used it: c*b_max/(ku*m_e c^2), ' // &
      'divided by sqrt(2) for a planar device. Zero off an FEL segment.', 'ele', aw, merr)
call fel_h5_real (g_id, 'l_period', 'm', 'period', 'Undulator period.', 'ele', l_per, merr)
call fel_h5_real (g_id, 'ku', '1/m', 'ku', 'Undulator wavenumber, twopi/l_period.', &
      'ele', ku, merr)
call fel_h5_flag (g_id, 'helical', 'helical', &
      'One for a helical device, zero for planar.', 'ele', helical, merr)
call fel_h5_real (g_id, 'k1', '1/m^2', 'k1', &
      'Quadrupole strength, signed as Bmad signs it.', 'ele', k1, merr)
call fel_h5_real (g_id, 'tilt', 'rad', 'tilt', &
      'Element tilt. On a planar FEL segment this is the wiggle-plane rotation, ' // &
      'which is the polarization spec.', 'ele', tilt, merr)
call fel_h5_real (g_id, 'z_offset', 'm', 'z offset', &
      'Longitudinal misalignment, the inter-segment phasing knob.', 'ele', z_off, merr)
call H5Gclose_f (g_id, h5e)

! The axis itself, into the coords/ group the stats writer left open for it. Every name
! in an @axes attribute resolves to a coords/ dataset, and lattice/ is the only table
! whose axis this routine knows the length of.

call H5Gopen_f (f_id, 'coords', g_id, h5e)
if (h5e < 0) then
  merr = .true.
else
  call fel_h5_int (g_id, 'ele', '1', 'element', &
        'Lattice element index, the axis the lattice/ table runs over. coords/ix_ele ' // &
        'is a variable whose values index it.', 'ele', [(ie, ie = 0, ne)], merr)
  call H5Gclose_f (g_id, h5e)
endif

if (merr) call out_io (s_warn$, r_name, 'Could not write all of lattice/.')

call H5Fclose_f (f_id, h5e)

end subroutine fel_write_lattice

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
