!+
! Program lucifer
!
! FEL tracker validated against Genesis 1.3 Version 4.
!
! The documentation is in docs/:
!   fel-physics.md       The physics: every equation the code integrates, its Genesis
!                         provenance, and the check that pins it.
!   input-reference.md   Every input key: default, meaning, and what refuses it.
!                         This is the normative reference. This header is a summary.
!   user-guide.md        Building, describing a run, running it, and the output files.
!   reading-output.md    What the output files hold and how to read them.
!   validation.md        What is checked against Genesis, how, and the recorded results.
!   performance.md       Where a run spends its time, on a named machine.
!
! The program reads one TOML input file through the schema (fel_schema_mod). Each table sets
! structs directly, as Tao's &tao_params does. Defaults are in the struct declarations.
!
!   lattice         The Bmad lattice file, the one top-level key.
!   [run]           The run switches, fel_global_struct.
!   [window]        The time window and the slice spacing, fel_slicing_struct.
!   [beam]          Bmad's beam_init bunch description, the load and its switches, the
!                   files, and [beam.resample], the resampler.
!   [field]         The radiation starting condition and the field files. The field
!                   record is the time window, so the window is set beside it.
!   [wake]          The chamber wake. [space_charge] the space-charge solver.
!   [bmad]          The bmad_com and space_charge_com switches the tracker honors.
!
! There are three ways in, all through the same input file:
!   A pair of Genesis dumps. Both codes then track from bitwise-identical initial conditions.
!   A quiet start generated from beam_init.
!   An imported particle distribution.
! Time dependence follows from the starting state alone. A multi-slice window makes a
! time-dependent run, and a single slice makes the steady state. There is no separate
! switch. A key the schema does not name is refused with the keys of its table, and a
! key in the wrong table with the table it belongs in. The arguments after the input file set
! keys of the input file, each table.key=value, read through the same schema after it.
!
! Usage:
!   lucifer <input.toml> [<table>.<key>=<value> ...]
!   lucifer --schema <schema.json>
!-

program lucifer

use fel_struct
use fel_input_mod
use fel_schema_mod
use fel_setup_mod
use fel_init_mod
use fel_track_line_mod
use fel_io_mod
use fel_timer_mod
use fel_device_mod, only: fel_device_query
use bmad_version_mod, only: bmad_version_date
use fel_version_mod, only: lucifer_version$
use, intrinsic :: iso_fortran_env, only: output_unit

implicit none

! The driver is read-parse-call (docs/user-guide.md). The input reader fills the input
! structs, and the library builds and walks the run. Every library error returns here.
! This is the one place that stops. The check instruments (split_weights, swap_beam_xy,
! gen_test_weights, resample_split_weights, wake.write_kernels) are carried in
! the input structs and act inside init and setup.

type (fel_run_struct), target :: run
integer n_arg, iu_k, i_k
logical err, dev_usable
character(400) param_file
character(400), allocatable :: overrides(:)
character(256) dev_reason
character(64) dev_name
character(*), parameter :: r_name = 'lucifer'

!

! The first argument is the input file. The examples name it lucifer.toml by convention. No name
! is looked for here. The arguments after it are overrides, table.key=value, which
! fel_read_input applies after the input file. --schema <file> writes the input file's schema as JSON
! Schema and stops. --version prints the Lucifer version and the Bmad version and stops.

! The run clock starts before anything is read. So the footer's run total covers the
! whole process, and the work outside the walk is the difference.

call fel_timer_reset()

n_arg = command_argument_count()
if (n_arg == 1) then
  call get_command_argument (1, param_file)
  if (param_file == '--version') then
    call out_io (s_blank$, r_name, 'Lucifer version ' // lucifer_version$)
    call out_io (s_blank$, r_name, 'Bmad version ' // bmad_version_date)
    flush (output_unit)
    stop
  endif
endif
if (n_arg == 2) then
  call get_command_argument (1, param_file)
  if (param_file == '--schema') then
    call get_command_argument (2, param_file)
    open (newunit = iu_k, file = trim(param_file), action = 'write', status = 'replace', iostat = i_k)
    if (i_k /= 0) then
      call out_io (s_error$, r_name, 'CANNOT WRITE THE SCHEMA TO: ' // trim(param_file))
      stop 1
    endif
    call fel_schema_write_json (iu_k)
    close (iu_k)
    call out_io (s_blank$, r_name, 'Schema version ' // fel_schema_version$ // ' written to ' // trim(param_file))
    stop
  endif
endif
if (n_arg == 0) then
  call out_io (s_blank$, r_name, [ &
        'Usage:                                                                  ', &
        '  lucifer <input.toml> [<table>.<key>=<value> ...]                      ', &
        '  lucifer --schema <schema.json>                                        ', &
        '  lucifer --version                                                     ', &
        '                                                                        ', &
        'One TOML input file names the lattice and holds the tables [run],       ', &
        '[window], [beam], [field], [wake], [space_charge] and [bmad]. Any file  ', &
        'name works. The examples name theirs lucifer.toml by convention, which  ', &
        'is a convention rather than a default looked for here. An argument after', &
        'the input file sets one key of it, beam.n_particle=4096 or              ', &
        'run.out_root=s7, read after the input file through the same schema.     ', &
        '--schema writes the schema as JSON Schema, and --version prints the two ', &
        'versions below. docs/ documents every key, every output file and the    ', &
        'lattice attributes.                                                     '])
  ! Separate from the block above. A character array constructor needs every element the
  ! same length, and the length of a version string is not known here.
  call out_io (s_blank$, r_name, 'Lucifer version ' // lucifer_version$)
  call out_io (s_blank$, r_name, 'Bmad version ' // bmad_version_date)
  ! Print what run.device can be set to in this build. The backend is a property of
  ! the build and the machine. Where there is none, the reason names the toolchain
  ! requirement and not the platform. The validation suite reads this line to decide whether
  ! the device checks can run.
  call fel_device_query (dev_usable, dev_name, dev_reason)
  if (dev_usable) then
    call out_io (s_blank$, r_name, 'Device backend: metal, on ' // trim(dev_name))
  else
    call out_io (s_blank$, r_name, 'Device backend: none. ' // trim(dev_reason))
  endif
  ! Flush before the stop. Otherwise the runtime's STOP line can reach the terminal first
  ! and the usage text appears to follow the exit.
  flush (output_unit)
  stop 1
endif
call get_command_argument (1, param_file)
allocate (overrides(n_arg-1))
do i_k = 2, n_arg
  call get_command_argument (i_k, overrides(i_k-1))
enddo

call fel_read_input (param_file, run, err, overrides)
if (err) stop 1

call fel_setup_lattice (run, err)
if (err) stop 1

call fel_init_beam (run, err)
if (err) stop 1

call fel_init_wavefront (run, err)
if (err) stop 1

call fel_setup_schedule (run, err)
if (err) stop 1

call fel_write_header (run)

! Check instrument: write the transcribed single-particle wake kernels for comparison
! with Genesis (fel-physics.md sec-wakes). fel_wake_init, called from
! fel_setup_schedule, builds the kernels.
! Note: The s = 0 entries carry the Bane self-slice half factor.

if (run%chamber_wake%write_kernels /= '' .and. run%coll%wake%on) then
  open (newunit = iu_k, file = trim(run%chamber_wake%write_kernels), action = 'write')
  write (iu_k, '(a)') '# s [m]   wakeres   wakegeo   wakerou   [eV/(m electron)]; s=0 rows are HALVED (Bane self-slice)'
  do i_k = 1, run%coll%wake%ns
    write (iu_k, '(4es24.15e3)') run%coll%wake%ds * (i_k-1), run%coll%wake%wakeres(i_k), &
                                 run%coll%wake%wakegeo(i_k), run%coll%wake%wakerou(i_k)
  enddo
  close (iu_k)
  call out_io (s_info$, r_name, 'Wrote wake kernels: ' // trim(run%chamber_wake%write_kernels))
endif

if (run%global%write_initial .or. run%global%load_only) then
  ! Both dumps use the same writers as the final ones. So the initial state is a file the
  ! tracker can read back.
  call fel_dump_field_set (run, trim(run%global%out_root) // '-initial', err, ix_done = run%i_start - 1)
  if (err) stop 1
  call fel_dump_beam (run, run%lat%branch(0)%ele(run%i_start), trim(run%global%out_root) // '-initial', err)
  if (err) stop 1
endif

! The pre-run. The header above has printed everything the run derives and everything
! it will cost. So the cost and the load of an input file can be read before any tracking is
! done (docs/user-guide.md).

if (run%global%load_only) then
  call out_io (s_info$, r_name, 'load_only is set. The header above is the pre-run: what the input file ' // &
        'derives, what it will cost and where its load stands. The initial state is written and ' // &
        'nothing is tracked.')
  stop 0
endif

call track_fel_line (run, err)
if (err) stop 1

! Final dumps. The field records are first unrotated to time order: position is of the
! time window holds record slice 1 + mod(is-1+first, nslice). Genesis's field writer
! does the same on the fly (fel-physics.md sec-slippage).

call fel_dump_beam (run, run%lat%branch(0)%ele(run%i_end), trim(run%global%out_root) // '-final', err, &
                    s_pos = run%lat%branch(0)%ele(run%i_end)%s)
if (err) stop 1

call fel_dump_field_set (run, trim(run%global%out_root) // '-final', err, ix_done = run%i_end)
if (err) stop 1

call fel_finalize_diagnostics (run, err)
if (err) stop 1

call fel_write_footer (run)

call wavefront_fft_free()

end program lucifer
