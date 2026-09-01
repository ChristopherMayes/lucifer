!+
! Program lucifer_smoke_test
!
! The library-contract check (doc/user-guide.md): drive the FEL tracker without the
! namelist layer. Fill the input structs in code, call the library, write the dumps.
! The harness compares its outputs dataset-identically against a namelist-driven run
! of the same configuration, and runs it twice in one process (both passes here) to
! prove re-entrancy: pass 2 must be bit-identical to pass 1 and to a separate process.
!
!   lucifer_smoke_test <lat_file> <out_root_pass1> <out_root_pass2>
!
! With an unreadable lattice file the library must return an error (never exit): this
! program then prints its proof-of-return line and stops with code 2. The harness's
! library-error check requires exactly that behavior.
!-

program lucifer_smoke_test

use fel_struct
use fel_setup_mod
use fel_init_mod
use fel_track_line_mod
use fel_io_mod

implicit none

type (fel_run_struct), target :: run
integer ipass, n_arg
logical err
character(400) lat_file, out_root(2)

!

n_arg = command_argument_count()
if (n_arg /= 3) then
  print '(a)', 'Usage: lucifer_smoke_test <lat_file> <out_root_pass1> <out_root_pass2>'
  stop 1
endif
call get_command_argument (1, lat_file)
call get_command_argument (2, out_root(1))
call get_command_argument (3, out_root(2))

do ipass = 1, 2

  ! A fresh run state per pass: the default constructor carries every declaration
  ! default, exactly as a freshly declared fel_run_struct.

  run = fel_run_struct()

  ! The configuration, set in code, no namelist anywhere. A small seeded
  ! steady-state case (the harness's twin namelist file carries the same values).

  run%lat_file = lat_file
  run%global%out_root = out_root(ipass)
  run%global%interlude_model = 'genesis'
  run%global%write_diag = .true.
  run%global%ran_seed = 777

  run%winit%lambda0 = 1e-10_rp
  run%winit%seed_power = 1e7_rp
  run%winit%seed_waist_size = 30e-6_rp
  run%winit%grid_n_pts = 63
  run%winit%grid_half_width = 2e-4_rp

  run%beam_init%n_particle = 2048
  run%beam_init%bunch_charge = 1.000692285594e-15_rp
  run%beam_init%sig_z = 0
  run%beam_init%sig_pz = 8.804506566858e-5_rp
  run%beam_init%a_norm_emit = 4e-7_rp
  run%beam_init%b_norm_emit = 4e-7_rp

  ! The library sequence. Every error returns. The proof-of-return line below is
  ! what the harness's library-error check greps for.

  err = .false.
  call fel_setup_lattice (run, err)
  if (.not. err) call fel_init_beam (run, err)
  if (.not. err) call fel_init_wavefront (run, err)
  if (.not. err) call fel_setup_schedule (run, err)
  if (.not. err) call track_fel_line (run, err)
  if (.not. err) call fel_dump_beam (run, run%lat%branch(0)%ele(run%i_end), trim(out_root(ipass)) // '-final', err)
  if (.not. err) call fel_dump_field_set (run, trim(out_root(ipass)) // '-final', err)
  if (.not. err) call fel_finalize_diagnostics (run, err)

  if (err) then
    print '(a)', 'lucifer_smoke_test: the library returned an error (no exit inside the library).'
    stop 2
  endif
  print '(a, i0, a)', 'lucifer_smoke_test: pass ', ipass, ' complete.'
enddo

call wavefront_fft_free()

end program lucifer_smoke_test
