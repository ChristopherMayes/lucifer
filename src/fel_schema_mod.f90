!+
! Module fel_schema_mod
!
! The schema of the input file: one table of every key an input file can set, with its TOML table,
! its type, its unit, its bounds and its description, and the structure component it
! fills. Everything that reads or describes the input file goes through this table. The reader
! (fel_input_mod) walks it, lucifer --schema writes it out as JSON Schema, the stats
! file's params/ group records every key from it, and the echo in meta/input_echo is
! the resolved input file written back from it. A key added here is read, refused, recorded,
! echoed and documented with no other edit.
!
! The defaults are not in this table. They are the declarations' own (fel_struct), and
! every default this module reports is read from a default-constructed run, so the table
! cannot disagree with the structures.
!
! The input file is one TOML file. Its tables are [run] for the run switches
! (fel_global_struct), [window] for the slicing, [beam] for the bunch and the load with
! [beam.resample] for the resampler, [field] for the radiation start, [wake] for the
! chamber wake, [space_charge] for the space-charge solver and [bmad] for the Bmad
! switches the tracker honors. The lattice is the one top-level key.
!-

module fel_schema_mod

use fel_struct
use tomlf, only: toml_table, toml_array, toml_key, get_value, len

implicit none

! The schema version, written into every stats file and the JSON Schema. It moves when a
! key is added, renamed, retyped or rebounded, so a reader can tell which input file a run read.

character(*), parameter :: fel_schema_version$ = '1'

! The kinds a key can have. A path is a string that names a file, written into params/
! only when set. A string is a short word, a mode or a name.

integer, parameter :: key_real$ = 1, key_int$ = 2, key_logical$ = 3, key_string$ = 4
integer, parameter :: key_path$ = 5, key_int_list$ = 6, key_string_list$ = 7

! The tables, in the order the echo and the JSON Schema write them.

character(16), parameter :: fel_schema_tables$(8) = [character(16) :: &
      'run', 'window', 'beam', 'beam.resample', 'field', 'wake', 'space_charge', 'bmad']

!+
! Struct fel_key_struct
!
! One key of the input file. lo and hi bound a real or integer key; a bound at huge is no
! bound, and lo_open or hi_open makes it exclusive. recorded is false for the one key
! that is not physics input, out_root, which every output file carries already.
!-

type fel_key_struct
  character(16) :: table = ''
  character(32) :: name = ''
  integer :: kind = 0
  character(16) :: unit = ''
  character(32) :: label = ''
  character(240) :: descrip = ''
  real(rp) :: lo = -huge(1.0_rp), hi = huge(1.0_rp)
  logical :: lo_open = .false., hi_open = .false.
  logical :: recorded = .true.
  logical :: required = .false.
end type

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_schema_table (keys)
!
! Routine to build the key table, in input file order. This is the one list of what an input file
! can say.
!
! Output:
!   keys  -- fel_key_struct(:), allocatable: Every key of the input file.
!-

subroutine fel_schema_table (keys)

type (fel_key_struct), allocatable :: keys(:)
integer n

!

allocate (keys(140))
n = 0

! The lattice.

call add ('', 'lattice', key_path$, '', 'lattice', &
      'The Bmad lattice file. The FEL segments are its wigglers with fel_method set.', required = .true.)

! [run]: fel_global_struct.

call add ('run', 'out_root', key_string$, '', 'output root', &
      'Root of every output file name.', recorded = .false.)
call add ('run', 'interlude_model', key_string$, '', 'interludes', &
      'Interlude transport: bmad (the seam) or genesis (transcribed).')
call add ('run', 'transport_model', key_string$, '', 'undulator transport', &
      'Transverse transport inside averaged FEL elements: bmad, Bmad''s own kernel, or genesis, ' // &
      'the transcribed maps for validation only.')
call add ('run', 'unaveraged_steps_per_period', key_int$, '', 'steps per period', &
      'Integration steps per undulator period in unaveraged elements. Below 10 the push does not converge.', lo = 10.0_rp)
call add ('run', 'unaveraged_ramp_periods', key_int$, '', 'ramp periods', &
      'Periods the field ramps over at each end of an unaveraged element. 0 keeps the default of 2, -1 is a hard edge.', lo = -1.0_rp)
call add ('run', 'source_model', key_string$, '', 'source', &
      'The FEL source model: deposit, the per-particle scatter, or coherent, the coherent-Gaussian source.')
call add ('run', 'source_filter', key_logical$, '', 'source filter', &
      'The angular filter on the source term, transcribed from Genesis4. Off to compare against a code without one.')
call add ('run', 'source_filter_angle', key_real$, 'rad', 'filter angle', &
      'The filter''s edge as an angle. 0 derives it from the beam and the gain.', lo = 0.0_rp)
call add ('run', 'source_filter_tolerance', key_real$, '', 'filter tolerance', &
      'The largest fraction of source intensity the filter may remove inside the protected angle.', &
      lo = 0.0_rp, hi = 1.0_rp, lo_open = .true., hi_open = .true.)
call add ('run', 'source_filter_xcut', key_real$, '', 'filter xcut', &
      'The edge in Genesis4''s own unit, the shifted grid index over ngrid. Validation only.', lo = 0.0_rp)
call add ('run', 'source_filter_ycut', key_real$, '', 'filter ycut', &
      'The edge in Genesis4''s own unit, the shifted grid index over ngrid. Validation only.', lo = 0.0_rp)
call add ('run', 'source_filter_width', key_real$, '', 'filter width', &
      'The sigmoid''s width as a fraction of its edge. Small is a sharp edge.', lo = 0.0_rp, lo_open = .true.)
call add ('run', 'source_filter_mutate', key_logical$, '', 'filter mutate', &
      'Self-test of the filter check: apply the sigmoid to the propagated field instead of the source.')
call add ('run', 'track_start', key_string$, '', 'track start', &
      'Element locator opening the walk. Blank is the first element.')
call add ('run', 'track_end', key_string$, '', 'track end', &
      'Element locator closing the walk. Blank is the last element.')
call add ('run', 'continuation', key_logical$, '', 'continuation', &
      'beam.file and field.file are a checkpoint pair the run continues from, not an initialization.')
call add ('run', 'comb_ds_save', key_real$, 'm', 'comb', &
      'Minimum z advance between per-record stats rows. 0 is every record, negative is none.')
call add ('run', 'interlude_ds_step', key_real$, 'm', 'interlude step', &
      'The length an interlude element is cut into, for stats rows and frames along it. 0 tracks it whole.')
call add ('run', 'dump_beam_at', key_string_list$, '', 'beam dumps', &
      'Element locators for mid-run beam dumps.')
call add ('run', 'dump_field_at', key_string_list$, '', 'field dumps', &
      'Element locators for mid-run field dumps.')
call add ('run', 'dump_at_comb', key_logical$, '', 'frame series', &
      'Write the beam and the field at every comb position as a frame series.')
call add ('run', 'dump_slice_first', key_int$, '', 'first frame slice', &
      'The first slice a frame carries, in window order.', lo = 1.0_rp)
call add ('run', 'dump_slice_last', key_int$, '', 'last frame slice', &
      'The last slice a frame carries. -1 is the last slice of the window.', lo = -1.0_rp)
call add ('run', 'dump_reduced', key_logical$, '', 'field reductions', &
      'Write the field''s own reductions per record into the stats file.')
call add ('run', 'ran_seed', key_int$, '', 'seed', &
      'The one random seed: generation, import and noise.')
call add ('run', 'write_diag', key_logical$, '', 'write diag', &
      'Write the Genesis-comparison text diag file, which is large.')
call add ('run', 'write_initial', key_logical$, '', 'write initial', &
      'Dump the initial state before tracking.')
call add ('run', 'load_only', key_logical$, '', 'load only', &
      'Build the initial state, dump it, stop.')
call add ('run', 'keep_escaped_field', key_logical$, '', 'keep escaped', &
      'Keep the light slippage carries out of the window, as the whole pulse at the exit. ' // &
      'On by default, since the light is what an FEL makes. A run whose window holds all ' // &
      'its light writes no file.')
call add ('run', 'write_escaped', key_logical$, '', 'write escaped', &
      'Also write the escaped file, each plane of the bank as it left the window. Off by ' // &
      'default, since the pulse holds the same light at the exit with each plane''s ' // &
      'transmission position and params. Needs keep_escaped_field.')
call add ('run', 'migrate', key_logical$, '', 'migrate', &
      'Move particles between slices when their phase leaves the slice window.')
call add ('run', 'migrate_check', key_logical$, '', 'migrate check', &
      'The migration''s bunching-invariance instrument.')
call add ('run', 'reference_run', key_logical$, '', 'reference run', &
      'No FEL interaction: Bmad tracks everything.')
call add ('run', 'record_environment', key_logical$, '', 'record environment', &
      'Record the login name and the working directory in meta/. Off, since a stats file travels.')
call add ('run', 'fp32_check', key_string$, '', 'fp32 check', &
      'The FP32 lockstep check: off, lockstep or freerun. For checks only.')
call add ('run', 'fp32_mutate', key_logical$, '', 'fp32 mutate', &
      'The FP32 check''s self-test: coarsen the residual.')
call add ('run', 'device', key_string$, '', 'device', &
      'The device backend: off, or metal on a build that carries it.')
call add ('run', 'device_timing', key_logical$, '', 'device timing', &
      'Per-pass device timing, one compute encoder a pass. Slows the run.')
call add ('run', 'device_dep_mutate', key_int$, '', 'deposit mutate', &
      'Shift the device deposit''s fixed-point scale, in bits. A check instrument.')

! [window]: fel_slicing_struct.

call add ('window', 'window_length', key_real$, 'm', 'window length', &
      'The time window. 0 derives it from the bunch and the slippage the line applies, or is one slice.', lo = 0.0_rp)
call add ('window', 'n_wavelength', key_int$, '', 'wavelengths per slice', &
      'Wavelengths per slice, Genesis4''s sample. 0 derives it from the gain.', lo = 0.0_rp)
call add ('window', 'n_slice', key_int$, '', 'slice count', &
      'The window as a slice count. 0 takes it from window_length. Both set is refused.', lo = 0.0_rp)
call add ('window', 'current', key_real$, 'A', 'flat current', &
      'A flat current. 0 derives the current from the bunch''s z structure.', lo = 0.0_rp)

! [beam]: beam_init_struct's honored components and the beam-side scalars.

call add ('beam', 'file', key_path$, '', 'beam file', &
      'An openPMD particle dump to start from, already in slices. Blank loads a bunch instead.')
call add ('beam', 'position_file', key_path$, '', 'position file', &
      'An openPMD-beamphysics bunch to load, read whole and then sliced by load_mode.')
call add ('beam', 'n_particle', key_int$, '', 'particles per slice', &
      'Macroparticles per slice after the phase copies, a positive multiple of beamlet_size. ' // &
      'In the keep mode, the particles of the generated bunch.', lo = 1.0_rp)
call add ('beam', 'bunch_charge', key_real$, 'C', 'charge', &
      'Total bunch charge. The current is derived from it and never input.', lo = 0.0_rp)
call add ('beam', 'sig_z', key_real$, 'm', 'bunch length', &
      'Bunch length, read with z_distribution. 0 is the steady state of one slice.', lo = 0.0_rp)
call add ('beam', 'sig_pz', key_real$, '', 'momentum spread', &
      'Fractional momentum spread dP/P0.', lo = 0.0_rp)
call add ('beam', 'a_norm_emit', key_real$, 'm rad', 'a norm emit', &
      'Normalized a-mode emittance.', lo = 0.0_rp)
call add ('beam', 'b_norm_emit', key_real$, 'm rad', 'b norm emit', &
      'Normalized b-mode emittance.', lo = 0.0_rp)
call add ('beam', 'z_distribution', key_string$, '', 'current profile', &
      'The current profile: ran_gauss, a Gaussian of sig_z, or grid, flat from z_min to z_max.')
call add ('beam', 'z_min', key_real$, 'm', 'z min', &
      'The flat profile''s lower z bound, read under grid only.')
call add ('beam', 'z_max', key_real$, 'm', 'z max', &
      'The flat profile''s upper z bound, read under grid only.')
call add ('beam', 'load_mode', key_string$, '', 'load mode', &
      'How the bunch becomes slices: sample, the same count in every slice with weights from the ' // &
      'slice charge, or keep, every particle kept with the counts following the charge.')
call add ('beam', 'quiet_start', key_logical$, '', 'quiet start', &
      'Load beamlets of beamlet_size phase copies, quiet below that harmonic.')
call add ('beam', 'beamlet_size', key_int$, '', 'beamlet size', &
      'Copies per beamlet, and the harmonics the noise resolves.', lo = 1.0_rp)
call add ('beam', 'shot_noise', key_logical$, '', 'shot noise', &
      'Impose the physical noise on the load''s groups. Time-dependent windows only.')
call add ('beam', 'write_genesis_dist', key_path$, '', 'genesis distribution', &
      'Write the bunch as a Genesis4 importdistribution file.')
call add ('beam', 'write_openpmd_file', key_path$, '', 'openPMD bunch', &
      'Write the bunch as openPMD-beamphysics.')
call add ('beam', 'split_weights', key_logical$, '', 'split weights', &
      'Check instrument: coincident copies carrying 1/3 and 2/3 of each weight after loading.')
call add ('beam', 'swap_beam_xy', key_logical$, '', 'swap xy', &
      'Check instrument: swap (x, px) with (y, py) after generation.')
call add ('beam', 'gen_test_weights', key_logical$, '', 'test weights', &
      'Check instrument: alternate beamlet weights of 0.25 and 1.75, charge preserving.')
call add ('beam', 'resample_split_weights', key_logical$, '', 'resample split', &
      'Check instrument: the split-weight copies applied to the bunch before it is sliced.')

! [beam.resample]: fel_resample_param_struct.

call add ('beam.resample', 'slice_width', key_real$, '', 'slice width', &
      'Sampling window over bunch length, Genesis4''s slicewidth.', lo = 0.0_rp, lo_open = .true.)
call add ('beam.resample', 'n_particle_per_slice', key_int$, '', 'particles', &
      'Macroparticles per slice on the validation route.', lo = 1.0_rp)
call add ('beam.resample', 'beamlet_size', key_int$, '', 'beamlet size', &
      'Beamlet size on the validation route.', lo = 1.0_rp)
call add ('beam.resample', 'n_slice', key_int$, '', 'slices', &
      'Slice count. 0 derives it from the bunch length and the spacing.', lo = 0.0_rp)
call add ('beam.resample', 'use_beam_init', key_logical$, '', 'use generated bunch', &
      'Validation route: resample a bunch Bmad generated from [beam] instead of the analytic loader.')

! [field]: wavefront_init_struct and the field files.

call add ('field', 'file', key_string_list$, '', 'field files', &
      'The openPMD wavefront file to start from, or one per harmonic member. Blank is a seed.')
call add ('field', 'lambda0', key_real$, 'm', 'wavelength', &
      'The resonant wavelength. Required for generation.', lo = 0.0_rp)
call add ('field', 'grid_n_pts', key_int$, '', 'grid points', &
      'Transverse grid points per side. 0 derives it from the beam.', lo = 0.0_rp)
call add ('field', 'grid_half_width', key_real$, 'm', 'grid half width', &
      'Transverse grid half width. 0 derives it from the beam.', lo = 0.0_rp)
call add ('field', 'seed_power', key_real$, 'W', 'seed power', &
      'Gaussian seed power. 0 is a dark start.', lo = 0.0_rp)
call add ('field', 'seed_waist_size', key_real$, 'm', 'seed waist', &
      'The seed''s intensity 1/e^2 radius.', lo = 0.0_rp)
call add ('field', 'seed_polarization', key_string$, '', 'seed polarization', &
      'The seed''s polarization, x or y.')
call add ('field', 'harmonics', key_int_list$, '', 'harmonics', &
      'The field set: 1 first, then further harmonic numbers in increasing order.')

! [wake]: fel_chamber_wake_init_struct.

call add ('wake', 'on', key_logical$, '', 'wake on', 'The chamber wake.')
call add ('wake', 'model', key_string$, '', 'wake model', &
      'Which chamber-wake implementation runs: genesis, the transcribed solver.')
call add ('wake', 'loss', key_real$, 'eV/m', 'external loss', 'External loss.')
call add ('wake', 'radius', key_real$, 'm', 'chamber radius', &
      'Chamber radius, or the half gap of a flat chamber.', lo = 0.0_rp, lo_open = .true.)
call add ('wake', 'conductivity', key_real$, '1/(Ohm m)', 'conductivity', &
      'DC conductivity. 0 is no resistive wake.', lo = 0.0_rp)
call add ('wake', 'relaxation', key_real$, 'm', 'relaxation', &
      'AC relaxation distance c tau.', lo = 0.0_rp)
call add ('wake', 'roundpipe', key_logical$, '', 'round pipe', &
      'A round chamber. False is flat, parallel plates.')
call add ('wake', 'material', key_string$, '', 'material', &
      'CU or AL, a shortcut for the conductivity and relaxation.')
call add ('wake', 'gap', key_real$, 'm', 'gap', 'Undulator gap. 0 is no geometric wake.', lo = 0.0_rp)
call add ('wake', 'lgap', key_real$, 'm', 'gap period', 'Period of the gaps.', lo = 0.0_rp, lo_open = .true.)
call add ('wake', 'hrough', key_real$, 'm', 'roughness', 'Roughness amplitude. 0 is no roughness wake.', lo = 0.0_rp)
call add ('wake', 'lrough', key_real$, 'm', 'roughness period', 'Roughness period.', lo = 0.0_rp, lo_open = .true.)
call add ('wake', 'write_kernels', key_path$, '', 'kernel file', &
      'Check instrument: write the transcribed kernels to this file.')

! [space_charge]: fel_space_charge_input_struct.

call add ('space_charge', 'model', key_string$, '', 'space-charge model', &
      'Which space-charge implementation runs: genesis, the transcribed solver.')
call add ('space_charge', 'rmax', key_real$, 'm', 'radial extent', &
      'Radial grid extent scale. It grows adaptively as Genesis4''s does.', lo = 0.0_rp)
call add ('space_charge', 'ngrid', key_int$, '', 'radial points', 'Radial grid points.', lo = 1.0_rp)
call add ('space_charge', 'nz', key_int$, '', 'harmonics', &
      'Longitudinal harmonics. 0 disables the short-range solve.', lo = 0.0_rp)
call add ('space_charge', 'nphi', key_int$, '', 'azimuthal modes', 'Azimuthal modes, -nphi to nphi.', lo = 0.0_rp)
call add ('space_charge', 'longrange', key_logical$, '', 'long range', 'The whole-window long-range term.')

! [bmad]: the Bmad switches the tracker honors, bmad_com and space_charge_com.

call add ('bmad', 'radiation_damping_on', key_logical$, '', 'radiation damping', 'Bmad''s radiation damping in the interludes.')
call add ('bmad', 'radiation_fluctuations_on', key_logical$, '', 'radiation fluctuations', 'Bmad''s radiation fluctuations in the interludes.')
call add ('bmad', 'csr_and_space_charge_on', key_logical$, '', 'csr and space charge', 'Bmad''s master switch for CSR and space charge.')
call add ('bmad', 'sr_wakes_on', key_logical$, '', 'short-range wakes', 'Bmad''s short-range element wakes.')
call add ('bmad', 'lr_wakes_on', key_logical$, '', 'long-range wakes', 'Bmad''s long-range element wakes.')
call add ('bmad', 'default_ds_step', key_real$, 'm', 'default step', 'Bmad''s default tracking step.', lo = 0.0_rp)
call add ('bmad', 'space_charge_mesh_size', key_int_list$, '', 'space-charge mesh', 'Bmad''s space-charge mesh, three counts.')
call add ('bmad', 'sc_min_in_bin', key_int$, '', 'minimum in bin', 'Bmad''s minimum particles per space-charge bin.', lo = 1.0_rp)
call add ('bmad', 'rel_tol_tracking', key_real$, '', 'tracking tolerance', 'Bmad''s relative tracking tolerance.', lo = 0.0_rp)
call add ('bmad', 'particle_sigma_cutoff', key_real$, '', 'sigma cutoff', 'Bmad''s particle sigma cutoff for the space-charge mesh.')

keys = keys(1:n)

!-------------------------------------------------------
contains

subroutine add (table, name, kind, unit, label, descrip, lo, hi, lo_open, hi_open, recorded, required)

character(*) table, name, unit, label, descrip
integer kind
real(rp), optional :: lo, hi
logical, optional :: lo_open, hi_open, recorded, required

!

n = n + 1
keys(n)%table = table
keys(n)%name = name
keys(n)%kind = kind
keys(n)%unit = unit
keys(n)%label = label
keys(n)%descrip = descrip
if (present(lo)) keys(n)%lo = lo
if (present(hi)) keys(n)%hi = hi
if (present(lo_open)) keys(n)%lo_open = lo_open
if (present(hi_open)) keys(n)%hi_open = hi_open
if (present(recorded)) keys(n)%recorded = recorded
if (present(required)) keys(n)%required = required

end subroutine add

end subroutine fel_schema_table

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_schema_set (run, key, tab, err_flag)
!
! Routine to read one key from the TOML table it lives in and store it in run. The
! type must be the key's own, and a real or integer key must lie inside its bounds.
! A refusal names the table, the key and what was expected.
!
! Input:
!   key       -- fel_key_struct: The key, present in tab.
!   tab       -- toml_table: The table holding it.
!
! Output:
!   run       -- fel_run_struct: The component the key fills is set.
!   err_flag  -- logical: Set on a wrong type, a value outside the bounds or a list
!                  longer than the component holds.
!-

subroutine fel_schema_set (run, key, tab, err_flag)

type (fel_run_struct), target :: run
type (fel_key_struct) key
type (toml_table), intent(inout) :: tab
type (toml_array), pointer :: arr
logical err_flag

real(rp) rv
integer iv, stat, i, n
logical lv
character(:), allocatable :: sv
character(400), allocatable :: slist(:)
integer, allocatable :: ilist(:)
character(60) where
character(*), parameter :: r_name = 'fel_schema_set'

!

err_flag = .true.
where = '[' // trim(key%table) // '] ' // trim(key%name)
if (key%table == '') where = trim(key%name)

select case (key%kind)
case (key_real$)
  call get_value (tab, trim(key%name), rv, stat = stat)
  if (stat /= 0) then
    call refuse ('expected a number')
    return
  endif
  if (outside(rv)) return
case (key_int$)
  call get_value (tab, trim(key%name), iv, stat = stat)
  if (stat /= 0) then
    call refuse ('expected an integer')
    return
  endif
  if (outside(real(iv, rp))) return
case (key_logical$)
  call get_value (tab, trim(key%name), lv, stat = stat)
  if (stat /= 0) then
    call refuse ('expected true or false')
    return
  endif
case (key_string$, key_path$)
  call get_value (tab, trim(key%name), sv, stat = stat)
  if (stat /= 0) then
    call refuse ('expected a quoted string')
    return
  endif
case (key_int_list$)
  ! One integer or a list of them, since one entry is a common case.
  call get_value (tab, trim(key%name), iv, stat = stat)
  if (stat == 0) then
    allocate (ilist(1))
    ilist(1) = iv
  else
    call get_value (tab, trim(key%name), arr, requested = .false.)
    if (.not. associated(arr)) then
      call refuse ('expected an integer or a list of them, [1, 3]')
      return
    endif
    allocate (ilist(len(arr)))
    do i = 1, len(arr)
      call get_value (arr, i, ilist(i), stat = stat)
      if (stat /= 0) then
        call refuse ('expected an integer or a list of them, [1, 3]')
        return
      endif
    enddo
  endif
case (key_string_list$)
  ! One string or a list of strings, since one file is the common case.
  call get_value (tab, trim(key%name), sv, stat = stat)
  if (stat == 0) then
    allocate (slist(1))
    slist(1) = sv
  else
    call get_value (tab, trim(key%name), arr, requested = .false.)
    if (.not. associated(arr)) then
      call refuse ('expected a quoted string or a list of them')
      return
    endif
    allocate (slist(len(arr)))
    do i = 1, len(arr)
      call get_value (arr, i, sv, stat = stat)
      if (stat /= 0) then
        call refuse ('expected a quoted string or a list of them')
        return
      endif
      slist(i) = sv
    enddo
  endif
end select

! The store. Each key names the component it fills.

select case (trim(key%table) // '.' // trim(key%name))
case ('.lattice');                         run%lat_file = sv
case ('run.out_root');                     run%global%out_root = sv
case ('run.interlude_model');              run%global%interlude_model = sv
case ('run.transport_model');              run%global%transport_model = sv
case ('run.unaveraged_steps_per_period');  run%global%unaveraged_steps_per_period = iv
case ('run.unaveraged_ramp_periods');      run%global%unaveraged_ramp_periods = iv
case ('run.source_model');                 run%global%source_model = sv
case ('run.source_filter');                run%global%source_filter = lv
case ('run.source_filter_angle');          run%global%source_filter_angle = rv
case ('run.source_filter_tolerance');      run%global%source_filter_tolerance = rv
case ('run.source_filter_xcut');           run%global%source_filter_xcut = rv
case ('run.source_filter_ycut');           run%global%source_filter_ycut = rv
case ('run.source_filter_width');          run%global%source_filter_width = rv
case ('run.source_filter_mutate');         run%global%source_filter_mutate = lv
case ('run.track_start');                  run%global%track_start = sv
case ('run.track_end');                    run%global%track_end = sv
case ('run.continuation');                 run%global%continuation = lv
case ('run.comb_ds_save');                 run%global%comb_ds_save = rv
case ('run.interlude_ds_step');            run%global%interlude_ds_step = rv
case ('run.dump_beam_at');                 if (fits(size(run%global%dump_beam_at))) run%global%dump_beam_at(1:size(slist)) = slist
case ('run.dump_field_at');                if (fits(size(run%global%dump_field_at))) run%global%dump_field_at(1:size(slist)) = slist
case ('run.dump_at_comb');                 run%global%dump_at_comb = lv
case ('run.dump_slice_first');             run%global%dump_slice_first = iv
case ('run.dump_slice_last');              run%global%dump_slice_last = iv
case ('run.dump_reduced');                 run%global%dump_reduced = lv
case ('run.ran_seed');                     run%global%ran_seed = iv
case ('run.write_diag');                   run%global%write_diag = lv
case ('run.write_initial');                run%global%write_initial = lv
case ('run.load_only');                    run%global%load_only = lv
case ('run.keep_escaped_field');           run%global%keep_escaped_field = lv;  run%bank_set = .true.
case ('run.write_escaped');                run%global%write_escaped = lv
case ('run.migrate');                      run%global%migrate = lv
case ('run.migrate_check');                run%global%migrate_check = lv
case ('run.reference_run');                run%global%reference_run = lv
case ('run.record_environment');           run%global%record_environment = lv
case ('run.fp32_check');                   run%global%fp32_check = sv
case ('run.fp32_mutate');                  run%global%fp32_mutate = lv
case ('run.device');                       run%global%device = sv
case ('run.device_timing');                run%global%device_timing = lv
case ('run.device_dep_mutate');            run%global%device_dep_mutate = iv
case ('window.window_length');             run%slicing%window_length = rv
case ('window.n_wavelength');              run%slicing%n_wavelength = iv
case ('window.n_slice');                   run%slicing%n_slice = iv
case ('window.current');                   run%slicing%current = rv
case ('beam.file');                        run%bparam%beam_file = sv
case ('beam.position_file');               run%beam_init%position_file = sv
case ('beam.n_particle');                  run%beam_init%n_particle = iv
case ('beam.bunch_charge');                run%beam_init%bunch_charge = rv
case ('beam.sig_z');                       run%beam_init%sig_z = rv
case ('beam.sig_pz');                      run%beam_init%sig_pz = rv
case ('beam.a_norm_emit');                 run%beam_init%a_norm_emit = rv
case ('beam.b_norm_emit');                 run%beam_init%b_norm_emit = rv
case ('beam.z_distribution');              run%beam_init%distribution_type(3) = sv
case ('beam.z_min');                       run%beam_init%grid(3)%x_min = rv
case ('beam.z_max');                       run%beam_init%grid(3)%x_max = rv
case ('beam.load_mode');                   run%bparam%load_mode = sv
case ('beam.quiet_start');                 run%bparam%quiet_start = lv
case ('beam.beamlet_size');                run%bparam%beamlet_size = iv
case ('beam.shot_noise');                  run%bparam%shot_noise = lv
case ('beam.write_genesis_dist');          run%bparam%write_genesis_dist = sv
case ('beam.write_openpmd_file');          run%bparam%write_openpmd_file = sv
case ('beam.split_weights');               run%bparam%split_weights = lv
case ('beam.swap_beam_xy');                run%bparam%swap_beam_xy = lv
case ('beam.gen_test_weights');            run%bparam%gen_test_weights = lv
case ('beam.resample_split_weights');      run%bparam%resample_split_weights = lv
case ('beam.resample.slice_width');        run%resample%slice_width = rv
case ('beam.resample.n_particle_per_slice'); run%resample%n_particle_per_slice = iv
case ('beam.resample.beamlet_size');       run%resample%beamlet_size = iv
case ('beam.resample.n_slice');            run%resample%n_slice = iv
case ('beam.resample.use_beam_init');      run%resample%use_beam_init = lv
case ('field.file');                       if (fits(size(run%field_file))) run%field_file(1:size(slist)) = slist
case ('field.lambda0');                    run%winit%lambda0 = rv
case ('field.grid_n_pts');                 run%winit%grid_n_pts = iv
case ('field.grid_half_width');            run%winit%grid_half_width = rv
case ('field.seed_power');                 run%winit%seed_power = rv
case ('field.seed_waist_size');            run%winit%seed_waist_size = rv
case ('field.seed_polarization');          run%winit%seed_polarization = sv
case ('field.harmonics')
  if (.not. fits(size(run%winit%harmonics))) return
  run%winit%harmonics = 0
  run%winit%harmonics(1:size(ilist)) = ilist
case ('wake.on');                          run%chamber_wake%on = lv
case ('wake.model');                       run%chamber_wake%model = sv
case ('wake.loss');                        run%chamber_wake%loss = rv
case ('wake.radius');                      run%chamber_wake%radius = rv
case ('wake.conductivity');                run%chamber_wake%conductivity = rv
case ('wake.relaxation');                  run%chamber_wake%relaxation = rv
case ('wake.roundpipe');                   run%chamber_wake%roundpipe = lv
case ('wake.material');                    run%chamber_wake%material = sv
case ('wake.gap');                         run%chamber_wake%gap = rv
case ('wake.lgap');                        run%chamber_wake%lgap = rv
case ('wake.hrough');                      run%chamber_wake%hrough = rv
case ('wake.lrough');                      run%chamber_wake%lrough = rv
case ('wake.write_kernels');               run%chamber_wake%write_kernels = sv
case ('space_charge.model');               run%space_charge%model = sv
case ('space_charge.rmax');                run%space_charge%rmax = rv
case ('space_charge.ngrid');               run%space_charge%ngrid = iv
case ('space_charge.nz');                  run%space_charge%nz = iv
case ('space_charge.nphi');                run%space_charge%nphi = iv
case ('space_charge.longrange');           run%space_charge%longrange = lv
case ('bmad.radiation_damping_on');        bmad_com%radiation_damping_on = lv
case ('bmad.radiation_fluctuations_on');   bmad_com%radiation_fluctuations_on = lv
case ('bmad.csr_and_space_charge_on');     bmad_com%csr_and_space_charge_on = lv
case ('bmad.sr_wakes_on');                 bmad_com%sr_wakes_on = lv
case ('bmad.lr_wakes_on');                 bmad_com%lr_wakes_on = lv
case ('bmad.default_ds_step');             bmad_com%default_ds_step = rv
case ('bmad.space_charge_mesh_size')
  if (size(ilist) /= 3) then
    call refuse ('expected three counts, [nx, ny, nz]')
    return
  endif
  space_charge_com%space_charge_mesh_size = ilist
case ('bmad.sc_min_in_bin');               space_charge_com%sc_min_in_bin = iv
case ('bmad.rel_tol_tracking');            space_charge_com%rel_tol_tracking = rv
case ('bmad.particle_sigma_cutoff');       space_charge_com%particle_sigma_cutoff = rv
case default
  call out_io (s_error$, r_name, 'THE SCHEMA NAMES A KEY THE STORE DOES NOT: ' // trim(where))
  return
end select

err_flag = .false.

!-------------------------------------------------------
contains

subroutine refuse (expected)
character(*) expected
call out_io (s_error$, r_name, 'IN THE INPUT FILE, ' // trim(where) // ': ' // expected // '.', trim(key%descrip))
end subroutine refuse

! outside reports a value past the key's bounds, and refuses it.

function outside (v) result (bad)
real(rp) v
logical bad
character(40) lo_text, hi_text
bad = .false.
if (key%lo > -huge(1.0_rp)) then
  bad = v < key%lo .or. (key%lo_open .and. v == key%lo)
endif
if (key%hi < huge(1.0_rp)) then
  bad = bad .or. v > key%hi .or. (key%hi_open .and. v == key%hi)
endif
if (.not. bad) return
lo_text = ''; hi_text = ''
if (key%lo > -huge(1.0_rp)) then
  lo_text = merge('above ', 'from  ', key%lo_open) // bound_text(key%lo)
endif
if (key%hi < huge(1.0_rp)) then
  hi_text = merge('below ', 'up to ', key%hi_open) // bound_text(key%hi)
endif
call out_io (s_error$, r_name, 'IN THE INPUT FILE, ' // trim(where) // ' = ' // bound_text(v) // &
      ' IS OUTSIDE ITS RANGE, WHICH IS ' // trim(lo_text) // ' ' // trim(hi_text) // '.', trim(key%descrip))
end function outside

! bound_text writes a bound as the input file would, an integer key's as an integer.

function bound_text (b) result (t)
real(rp) b
character(:), allocatable :: t
if (key%kind == key_int$) then
  t = int_str(nint(b))
else
  t = real_text(b)
endif
end function bound_text

! fits reports whether a list fits the component, and refuses one that does not.

function fits (n_max) result (ok)
integer n_max, n_got
logical ok
if (allocated(slist)) then
  n_got = size(slist)
else
  n_got = size(ilist)
endif
ok = n_got <= n_max
if (.not. ok) then
  call out_io (s_error$, r_name, 'IN THE INPUT FILE, ' // trim(where) // ' HOLDS ' // int_str(n_got) // &
        ' ENTRIES AND AT MOST ' // int_str(n_max) // ' ARE TAKEN.')
endif
end function fits

end subroutine fel_schema_set

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_schema_value (run, key, rv, iv, lv, sv, ilist, slist)
!
! Routine to fetch one key's value from run into the argument of its kind. This is the
! one place that knows which component each key fills on the way out, as
! fel_schema_set is on the way in. A string list comes back without its blank tail.
!
! Input:
!   run    -- fel_run_struct: The run state.
!   key    -- fel_key_struct: The key.
!
! Output:
!   rv     -- real(rp): The value of a real key.
!   iv     -- integer: The value of an integer key.
!   lv     -- logical: The value of a logical key.
!   sv     -- character(:), allocatable: The value of a string or path key.
!   ilist  -- integer(:), allocatable: The value of an integer-list key.
!   slist  -- character(400), allocatable: The value of a string-list key.
!-

subroutine fel_schema_value (run, key, rv, iv, lv, sv, ilist, slist)

type (fel_run_struct), target :: run
type (fel_key_struct) key
real(rp) rv
integer iv, n
logical lv
character(:), allocatable :: sv
integer, allocatable :: ilist(:)
character(400), allocatable :: slist(:)

!

rv = 0; iv = 0; lv = .false.; sv = ''

select case (trim(key%table) // '.' // trim(key%name))
case ('.lattice');                         sv = run%lat_file
case ('run.out_root');                     sv = run%global%out_root
case ('run.interlude_model');              sv = run%global%interlude_model
case ('run.transport_model');              sv = run%global%transport_model
case ('run.unaveraged_steps_per_period');  iv = run%global%unaveraged_steps_per_period
case ('run.unaveraged_ramp_periods');      iv = run%global%unaveraged_ramp_periods
case ('run.source_model');                 sv = run%global%source_model
case ('run.source_filter');                lv = run%global%source_filter
case ('run.source_filter_angle');          rv = run%global%source_filter_angle
case ('run.source_filter_tolerance');      rv = run%global%source_filter_tolerance
case ('run.source_filter_xcut');           rv = run%global%source_filter_xcut
case ('run.source_filter_ycut');           rv = run%global%source_filter_ycut
case ('run.source_filter_width');          rv = run%global%source_filter_width
case ('run.source_filter_mutate');         lv = run%global%source_filter_mutate
case ('run.track_start');                  sv = run%global%track_start
case ('run.track_end');                    sv = run%global%track_end
case ('run.continuation');                 lv = run%global%continuation
case ('run.comb_ds_save');                 rv = run%global%comb_ds_save
case ('run.interlude_ds_step');            rv = run%global%interlude_ds_step
case ('run.dump_beam_at');                 slist = run%global%dump_beam_at
case ('run.dump_field_at');                slist = run%global%dump_field_at
case ('run.dump_at_comb');                 lv = run%global%dump_at_comb
case ('run.dump_slice_first');             iv = run%global%dump_slice_first
case ('run.dump_slice_last');              iv = run%global%dump_slice_last
case ('run.dump_reduced');                 lv = run%global%dump_reduced
case ('run.ran_seed');                     iv = run%global%ran_seed
case ('run.write_diag');                   lv = run%global%write_diag
case ('run.write_initial');                lv = run%global%write_initial
case ('run.load_only');                    lv = run%global%load_only
case ('run.keep_escaped_field');           lv = run%global%keep_escaped_field
case ('run.write_escaped');                lv = run%global%write_escaped
case ('run.migrate');                      lv = run%global%migrate
case ('run.migrate_check');                lv = run%global%migrate_check
case ('run.reference_run');                lv = run%global%reference_run
case ('run.record_environment');           lv = run%global%record_environment
case ('run.fp32_check');                   sv = run%global%fp32_check
case ('run.fp32_mutate');                  lv = run%global%fp32_mutate
case ('run.device');                       sv = run%global%device
case ('run.device_timing');                lv = run%global%device_timing
case ('run.device_dep_mutate');            iv = run%global%device_dep_mutate
case ('window.window_length');             rv = run%slicing%window_length
case ('window.n_wavelength');              iv = run%slicing%n_wavelength
case ('window.n_slice');                   iv = run%slicing%n_slice
case ('window.current');                   rv = run%slicing%current
case ('beam.file');                        sv = run%bparam%beam_file
case ('beam.position_file');               sv = run%beam_init%position_file
case ('beam.n_particle');                  iv = run%beam_init%n_particle
case ('beam.bunch_charge');                rv = run%beam_init%bunch_charge
case ('beam.sig_z');                       rv = run%beam_init%sig_z
case ('beam.sig_pz');                      rv = run%beam_init%sig_pz
case ('beam.a_norm_emit');                 rv = run%beam_init%a_norm_emit
case ('beam.b_norm_emit');                 rv = run%beam_init%b_norm_emit
case ('beam.z_distribution');              sv = run%beam_init%distribution_type(3)
case ('beam.z_min');                       rv = run%beam_init%grid(3)%x_min
case ('beam.z_max');                       rv = run%beam_init%grid(3)%x_max
case ('beam.load_mode');                   sv = run%bparam%load_mode
case ('beam.quiet_start');                 lv = run%bparam%quiet_start
case ('beam.beamlet_size');                iv = run%bparam%beamlet_size
case ('beam.shot_noise');                  lv = run%bparam%shot_noise
case ('beam.write_genesis_dist');          sv = run%bparam%write_genesis_dist
case ('beam.write_openpmd_file');          sv = run%bparam%write_openpmd_file
case ('beam.split_weights');               lv = run%bparam%split_weights
case ('beam.swap_beam_xy');                lv = run%bparam%swap_beam_xy
case ('beam.gen_test_weights');            lv = run%bparam%gen_test_weights
case ('beam.resample_split_weights');      lv = run%bparam%resample_split_weights
case ('beam.resample.slice_width');        rv = run%resample%slice_width
case ('beam.resample.n_particle_per_slice'); iv = run%resample%n_particle_per_slice
case ('beam.resample.beamlet_size');       iv = run%resample%beamlet_size
case ('beam.resample.n_slice');            iv = run%resample%n_slice
case ('beam.resample.use_beam_init');      lv = run%resample%use_beam_init
case ('field.file');                       slist = run%field_file
case ('field.lambda0');                    rv = run%winit%lambda0
case ('field.grid_n_pts');                 iv = run%winit%grid_n_pts
case ('field.grid_half_width');            rv = run%winit%grid_half_width
case ('field.seed_power');                 rv = run%winit%seed_power
case ('field.seed_waist_size');            rv = run%winit%seed_waist_size
case ('field.seed_polarization');          sv = run%winit%seed_polarization
case ('field.harmonics')
  n = count(run%winit%harmonics > 0)
  ilist = run%winit%harmonics(1:n)
case ('wake.on');                          lv = run%chamber_wake%on
case ('wake.model');                       sv = run%chamber_wake%model
case ('wake.loss');                        rv = run%chamber_wake%loss
case ('wake.radius');                      rv = run%chamber_wake%radius
case ('wake.conductivity');                rv = run%chamber_wake%conductivity
case ('wake.relaxation');                  rv = run%chamber_wake%relaxation
case ('wake.roundpipe');                   lv = run%chamber_wake%roundpipe
case ('wake.material');                    sv = run%chamber_wake%material
case ('wake.gap');                         rv = run%chamber_wake%gap
case ('wake.lgap');                        rv = run%chamber_wake%lgap
case ('wake.hrough');                      rv = run%chamber_wake%hrough
case ('wake.lrough');                      rv = run%chamber_wake%lrough
case ('wake.write_kernels');               sv = run%chamber_wake%write_kernels
case ('space_charge.model');               sv = run%space_charge%model
case ('space_charge.rmax');                rv = run%space_charge%rmax
case ('space_charge.ngrid');               iv = run%space_charge%ngrid
case ('space_charge.nz');                  iv = run%space_charge%nz
case ('space_charge.nphi');                iv = run%space_charge%nphi
case ('space_charge.longrange');           lv = run%space_charge%longrange
case ('bmad.radiation_damping_on');        lv = bmad_com%radiation_damping_on
case ('bmad.radiation_fluctuations_on');   lv = bmad_com%radiation_fluctuations_on
case ('bmad.csr_and_space_charge_on');     lv = bmad_com%csr_and_space_charge_on
case ('bmad.sr_wakes_on');                 lv = bmad_com%sr_wakes_on
case ('bmad.lr_wakes_on');                 lv = bmad_com%lr_wakes_on
case ('bmad.default_ds_step');             rv = bmad_com%default_ds_step
case ('bmad.space_charge_mesh_size');      ilist = space_charge_com%space_charge_mesh_size
case ('bmad.sc_min_in_bin');               iv = space_charge_com%sc_min_in_bin
case ('bmad.rel_tol_tracking');            rv = space_charge_com%rel_tol_tracking
case ('bmad.particle_sigma_cutoff');       rv = space_charge_com%particle_sigma_cutoff
case default
  call out_io (s_error$, 'fel_schema_value', 'THE SCHEMA NAMES A KEY THE VALUE GETTER DOES NOT: ' // &
        trim(key%table) // '.' // trim(key%name))
end select

! A string list is cut at its last set entry.

if (allocated(slist)) then
  n = 0
  do iv = 1, size(slist)
    if (slist(iv) /= '') n = iv
  enddo
  slist = slist(1:n)
  iv = 0
endif

end subroutine fel_schema_value

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_schema_text (run, key) result (txt)
!
! Routine to render one key's value in run as TOML text, the form an input file writes it in.
! The echo and the JSON Schema's defaults are built from this.
!
! Input:
!   run  -- fel_run_struct: The run state.
!   key  -- fel_key_struct: The key.
!
! Output:
!   txt  -- character(:), allocatable: The value as TOML text.
!-

function fel_schema_text (run, key) result (txt)

type (fel_run_struct), target :: run
type (fel_key_struct) key
character(:), allocatable :: txt
real(rp) rv
integer iv, j
logical lv
character(:), allocatable :: sv
integer, allocatable :: ilist(:)
character(400), allocatable :: slist(:)

!

call fel_schema_value (run, key, rv, iv, lv, sv, ilist, slist)

select case (key%kind)
case (key_real$);     txt = real_text(rv)
case (key_int$);      txt = int_str(iv)
case (key_logical$)
  if (lv) then
    txt = 'true'
  else
    txt = 'false'
  endif
case (key_string$, key_path$);  txt = fel_schema_quoted(sv)
case (key_int_list$)
  txt = '['
  do j = 1, size(ilist)
    if (j > 1) txt = txt // ', '
    txt = txt // int_str(ilist(j))
  enddo
  txt = txt // ']'
case (key_string_list$)
  txt = '['
  do j = 1, size(slist)
    if (j > 1) txt = txt // ', '
    txt = txt // fel_schema_quoted(slist(j))
  enddo
  txt = txt // ']'
end select

end function fel_schema_text

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_schema_quoted (s) result (q)
!
! Routine to quote a string as a TOML basic string, escaping the quote and the backslash.
!
! Input:
!   s  -- character(*): The string.
!
! Output:
!   q  -- character(:), allocatable: The quoted string.
!-

function fel_schema_quoted (s) result (q)

character(*) s
character(:), allocatable :: q
integer j

!

q = '"'
do j = 1, len_trim(s)
  select case (s(j:j))
  case ('"');   q = q // '\"'
  case ('\');  q = q // '\\'
  case default; q = q // s(j:j)
  end select
enddo
q = q // '"'

end function fel_schema_quoted

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function real_text (v) result (txt)
!
! Routine to render a real as TOML text: an integer value as "1.0", else in
! exponent form with the fewest significant digits, 15 to 17, that read back to the
! same value and the trailing zeros of the mantissa removed, so 1e-10 reads back as
! 1.0e-10, 0.75 as 7.5e-1 and 0.05 as 5.0e-2.
!
! Input:
!   v    -- real(rp): The value.
!
! Output:
!   txt  -- character(:), allocatable: The value as TOML text.
!-

function real_text (v) result (txt)

real(rp) v, back
character(:), allocatable :: txt
character(32) buf
character(12) fmt
integer ie, j, nd

!

if (v == 0) then
  txt = '0.0'
  return
endif

! The fewest significant digits, 15 to 17, that read back to the same value, so 0.05
! reads 5.0e-2 and not 5.0000000000000003e-2, while a value that needs all 17 keeps
! them.

do nd = 15, 17
  write (fmt, '(a, i0, a)') '(es24.', nd - 1, 'e3)'
  write (buf, fmt) v
  read (buf, *) back
  if (back == v) exit
enddo
buf = adjustl(buf)
ie = index(buf, 'E')
! Strip trailing zeros of the mantissa, keeping one digit after the point, and the
! leading zeros of the exponent, so 0.75 reads 7.5e-1.
j = ie - 1
do while (j > 1 .and. buf(j:j) == '0' .and. buf(j-1:j-1) /= '.')
  j = j - 1
enddo
txt = buf(1:j) // 'e' // exponent_text(buf(ie+1:))

!-------------------------------------------------------
contains

function exponent_text (e) result (t)
character(*) e
character(:), allocatable :: t
integer k
k = 2
do while (k < len_trim(e) .and. e(k:k) == '0')
  k = k + 1
enddo
if (e(1:1) == '-') then
  t = '-' // trim(e(k:))
else
  t = trim(e(k:))
endif
end function exponent_text

end function real_text

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_schema_write_toml (run, iu)
!
! Routine to write the resolved input file of run as TOML to an open unit, every key of the
! schema with its value, defaults made explicit. This is the echo in the stats file's
! meta/input_echo, and what an embedding program writes to save a configuration built
! in code.
!
! Input:
!   run -- fel_run_struct: The run state after reading.
!   iu  -- integer: Open unit to write to.
!-

subroutine fel_schema_write_toml (run, iu)

type (fel_run_struct), target :: run
type (fel_key_struct), allocatable :: keys(:)
integer iu, i, it
character(16) table

!

call fel_schema_table (keys)

do i = 1, size(keys)
  if (keys(i)%table /= '') cycle
  write (iu, '(3a)') trim(keys(i)%name), ' = ', fel_schema_text(run, keys(i))
enddo

do it = 1, size(fel_schema_tables$)
  table = fel_schema_tables$(it)
  write (iu, '(a)') ''
  write (iu, '(3a)') '[', trim(table), ']'
  do i = 1, size(keys)
    if (keys(i)%table /= table) cycle
    write (iu, '(3a)') trim(keys(i)%name), ' = ', fel_schema_text(run, keys(i))
  enddo
enddo

end subroutine fel_schema_write_toml

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_schema_write_json (iu)
!
! Routine to write the schema as JSON Schema, draft 2020-12, to an open unit. Each table
! is an object that allows no key the schema does not name, and each key carries its
! type, its default as the structures declare it, its description with its unit, and
! its bounds. An editor that reads JSON Schema validates an input file against it as it is
! typed, and the input reference's key tables are generated from it.
!
! Input:
!   iu  -- integer: Open unit to write to.
!-

subroutine fel_schema_write_json (iu)

type (fel_run_struct), target :: def
type (fel_key_struct), allocatable :: keys(:)
integer iu, i, it, n_in
logical key_kind_int
character(16) table
character(:), allocatable :: descrip

!

call fel_schema_table (keys)

write (iu, '(a)') '{'
write (iu, '(a)') '  "$schema": "https://json-schema.org/draft/2020-12/schema",'
write (iu, '(a)') '  "title": "Lucifer input",'
write (iu, '(a)') '  "description": "One TOML file: the lattice and the tables run, window, beam, field, wake, space_charge and bmad.",'
write (iu, '(3a)') '  "x-schema-version": "', fel_schema_version$, '",'
write (iu, '(a)') '  "type": "object",'
write (iu, '(a)') '  "additionalProperties": false,'
write (iu, '(a)') '  "required": ["lattice"],'
write (iu, '(a)') '  "properties": {'

do i = 1, size(keys)
  if (keys(i)%table /= '') cycle
  call key_json (keys(i), '    ', .true.)
enddo

do it = 1, size(fel_schema_tables$)
  table = fel_schema_tables$(it)
  if (index(table, '.') > 0) cycle    ! A sub-table is written inside its parent.
  write (iu, '(3a)') '    "', trim(table), '": {'
  write (iu, '(a)') '      "type": "object",'
  write (iu, '(a)') '      "additionalProperties": false,'
  write (iu, '(a)') '      "properties": {'
  n_in = count(keys%table == table) + count(keys%table == trim(table) // '.resample')
  do i = 1, size(keys)
    if (keys(i)%table /= table) cycle
    n_in = n_in - 1
    call key_json (keys(i), '        ', n_in > 0)
  enddo
  if (table == 'beam') then
    write (iu, '(a)') '        "resample": {'
    write (iu, '(a)') '          "type": "object",'
    write (iu, '(a)') '          "additionalProperties": false,'
    write (iu, '(a)') '          "properties": {'
    n_in = count(keys%table == 'beam.resample')
    do i = 1, size(keys)
      if (keys(i)%table /= 'beam.resample') cycle
      n_in = n_in - 1
      call key_json (keys(i), '            ', n_in > 0)
    enddo
    write (iu, '(a)') '          }'
    write (iu, '(a)') '        }'
  endif
  write (iu, '(a)') '      }'
  if (it < size(fel_schema_tables$)) then
    write (iu, '(a)') '    },'
  else
    write (iu, '(a)') '    }'
  endif
enddo

write (iu, '(a)') '  }'
write (iu, '(a)') '}'

!-------------------------------------------------------
contains

subroutine key_json (key, indent, comma)

type (fel_key_struct) key
character(*) indent
logical comma
character(:), allocatable :: jtype, items, bounds, dflt

!

key_kind_int = key%kind == key_int$
select case (key%kind)
case (key_real$);         jtype = 'number'
case (key_int$);          jtype = 'integer'
case (key_logical$);      jtype = 'boolean'
case (key_string$, key_path$); jtype = 'string'
case (key_int_list$);     jtype = 'array'; items = ', "items": {"type": "integer"}'
case (key_string_list$);  jtype = 'array'; items = ', "items": {"type": "string"}'
end select
if (.not. allocated(items)) items = ''

bounds = ''
if (key%lo > -huge(1.0_rp)) then
  if (key%lo_open) then
    bounds = bounds // ', "exclusiveMinimum": ' // bound_json(key%lo)
  else
    bounds = bounds // ', "minimum": ' // bound_json(key%lo)
  endif
endif
if (key%hi < huge(1.0_rp)) then
  if (key%hi_open) then
    bounds = bounds // ', "exclusiveMaximum": ' // bound_json(key%hi)
  else
    bounds = bounds // ', "maximum": ' // bound_json(key%hi)
  endif
endif

descrip = trim(key%descrip)
if (key%unit /= '') descrip = descrip // ' [' // trim(key%unit) // ']'
dflt = fel_schema_text(def, key)
! A default the input file would write as a list of one string stays a list in JSON too.
write (iu, '(a)', advance = 'no') indent // '"' // trim(key%name) // '": {"type": "' // jtype // '"' // items // &
      ', "default": ' // dflt // bounds // ', "description": "' // descrip // '"}'
if (comma) then
  write (iu, '(a)') ','
else
  write (iu, '(a)') ''
endif

end subroutine key_json

! bound_json writes a bound as a JSON number, an integer key's as an integer.

function bound_json (b) result (t)
real(rp) b
character(:), allocatable :: t
if (key_kind_int) then
  t = int_str(nint(b))
else
  t = real_text(b)
endif
end function bound_json

end subroutine fel_schema_write_json

end module fel_schema_mod
