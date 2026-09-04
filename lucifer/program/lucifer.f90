!+
! Program lucifer
!
! FEL tracker validated against Genesis 1.3 Version 4.
!
! The documentation is in lucifer/doc/:
!   fel-physics.md      The physics: every equation the code integrates, each piece's
!                        Genesis provenance, and which check pins it.
!   input-reference.md   Every namelist parameter: default, meaning, and what refuses it.
!                        The normative reference. This header is a summary only.
!   user-guide.md        Building, describing a run, running it, the output inventory.
!   reading-output.md    What the output files hold and how to read them.
!   validation.md        The keystone rule, the tier table, and the measured levels.
!   performance.md       Where a run spends its time, measured, with the machine named.
!
! The program reads three namelist groups from one input file, each setting structs
! directly (Tao's &tao_params pattern). Defaults live in the struct declarations.
!
!   &fel_params           lat_file, the global%... run switches, Bmad's own bmad_com and
!                         space_charge_com, and the chamber_wake%/space_charge% collective descriptions.
!   &fel_beam_init        Bmad's beam_init%... bunch description, the resample%... resampler,
!                         source and output files, and the beam-side check knobs.
!   &fel_wavefront_init   wavefront_init%... radiation starting condition (the field
!                         record is the time window, so the window lives here) and the
!                         field_file imports.
!
! Three ways in, all through the same input file: a pair of Genesis dumps (both codes
! then track from bitwise-identical initial conditions), a self-generated quiet start
! from beam_init, or an imported particle distribution. Time dependence follows from the
! starting state alone: a multi-slice window makes a time-dependent run, a single slice
! the steady state, no separate switch. The retired flat &fel_track_params group is
! refused, each parameter mapped to its new home.
!
! Usage:
!   lucifer <input_file>
!-

program lucifer

use fel_struct
use fel_input_mod
use fel_setup_mod
use fel_init_mod
use fel_track_line_mod
use fel_io_mod
use fel_timer_mod
use fel_device_mod, only: fel_device_query
use bmad_version_mod, only: bmad_version_date
use, intrinsic :: iso_fortran_env, only: output_unit

implicit none

! The driver is read-parse-call (doc/user-guide.md): the namelist layer fills the
! input structs, the library builds and walks the run, and every library error
! Returns here. This is the one place that stops. The check instruments
! (split_weights, swap_beam_xy, gen_test_weights, resample_split_weights,
! chamber_wake%write_kernels) ride along in the input structs and act inside init/setup.

type (fel_run_struct), target :: run
integer n_arg, iu_k, i_k
logical err, dev_usable
character(400) param_file
character(256) dev_reason
character(64) dev_name
character(*), parameter :: r_name = 'lucifer'

!

! The one argument is the input file. Its name is the caller's choice, and lucifer.in is
! the convention the examples follow rather than a name looked for here.

! The run clock starts before anything is read, so the footer's run total covers the
! whole process and the work outside the walk reads off as the difference.

call fel_timer_reset()

n_arg = command_argument_count()
if (n_arg /= 1) then
  call out_io (s_blank$, r_name, [ &
        'Usage:                                                                      ', &
        '  lucifer <input_file>                                                      ', &
        '                                                                            ', &
        'One input file names the lattice and holds the three namelist groups:       ', &
        '&fel_params, &fel_beam_init and &fel_wavefront_init. Any file name works.   ', &
        'The examples name theirs lucifer.in by convention, which is a convention    ', &
        'rather than a default looked for here. lucifer/doc/ documents every         ', &
        'parameter, every output file and the lattice attributes.                    '])
  ! Separate from the block above because a character array constructor needs every
  ! element the same length, and the version string's length is not this file's to know.
  call out_io (s_blank$, r_name, 'Bmad version ' // bmad_version_date)
  ! What global%device can be set to in this build, which otherwise takes a run to
  ! discover: the backend is a build-time and machine-time property, and the reason a
  ! build has none names the toolchain requirement rather than the platform. A test
  ! harness reads this line to decide whether the device checks can run at all.
  call fel_device_query (dev_usable, dev_name, dev_reason)
  if (dev_usable) then
    call out_io (s_blank$, r_name, 'Device backend: metal, on ' // trim(dev_name))
  else
    call out_io (s_blank$, r_name, 'Device backend: none. ' // trim(dev_reason))
  endif
  ! Flushed before the stop, or the runtime's own STOP line reaches the terminal first
  ! and the usage text appears to follow the exit.
  flush (output_unit)
  stop 1
endif
call get_command_argument (1, param_file)

call fel_read_input (param_file, run, err)
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

! Check instrument: export the transcribed single-particle wake kernels for
! cross-validation against Genesis (fel-physics.md sec-wakes). The kernels are built by
! fel_setup_schedule's fel_wake_init. Note: The s = 0 entries carry the Bane
! self-slice half factor.

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
  ! Both dumps go through the same writers as the final ones, so the initial state is a
  ! file the tracker can read straight back.
  call fel_dump_field_set (run, trim(run%global%out_root) // '-initial', err)
  if (err) stop 1
  call fel_dump_beam (run, run%lat%branch(0)%ele(run%i_start), trim(run%global%out_root) // '-initial', err)
  if (err) stop 1
endif

if (run%global%load_only) then
  call out_io (s_info$, r_name, 'load_only set; initial state written, no tracking.')
  stop 0
endif

call track_fel_line (run, err)
if (err) stop 1

! Final dumps. The field records are unrotated to time order first (position is of
! the time window holds record slice 1 + mod(is-1+first, nslice)), which is what
! Genesis's field writer does on the fly (fel-physics.md sec-slippage).

call fel_dump_beam (run, run%lat%branch(0)%ele(run%i_end), trim(run%global%out_root) // '-final', err)
if (err) stop 1

call fel_dump_field_set (run, trim(run%global%out_root) // '-final', err)
if (err) stop 1

call fel_finalize_diagnostics (run, err)
if (err) stop 1

call fel_write_footer (run)

call wavefront_fft_free()

end program lucifer
