!+
! Module fel_input_mod
!
! The namelist layer of the FEL tracker, kept apart the way Tao keeps its tao_init_*
! files. There are three groups, all read from one input file. Each group fills the
! structs of fel_struct directly, as Tao's &tao_params does. EG: global%out_root,
! bmad_com%radiation_damping_on, chamber_wake%radius, space_charge%nz, beam_init%n_particle
! and wavefront_init%lambda0 are set by component. This module is library code. An
! embedding program may use the parsing, or skip it and fill the structs in code.
!
!   &fel_params          lat_file, global, bmad_com, space_charge_com, wake, sc,
!                        chamber_wake%write_kernels
!   &fel_beam_init       beam_init, resample, beam_file, write_genesis_dist,
!                        write_openpmd_file, load_mode, quiet_start, beamlet_size, shot_noise,
!                        split_weights, swap_beam_xy, gen_test_weights,
!                        resample_split_weights
!   &fel_wavefront_init  wavefront_init, field_file
!
! A group that is absent keeps its defaults. A group that is present must parse.
! The retired flat &fel_track_params group is refused, and the error names the three
! groups that replace it.
!-

module fel_input_mod

use fel_struct

implicit none

contains

!------------------------------------------------------------------------------
!+
! Subroutine fel_read_input (param_file, run, err_flag)
!
! Routine to read the three namelist groups of param_file into the input structs of run.
! The namelist sets bmad_com and space_charge_com directly, since they are Bmad's own
! globals.
!
! Input:
!   param_file   -- character(*): The input file.
!
! Output:
!   run          -- fel_run_struct: Input structs filled (global, winit, wake_init,
!                     sc_init, beam_init, imp, bparam, lat_file, field_file). Nothing else
!                     is touched.
!   err_flag     -- logical: Set on any parse error or on a retired input name.
!-

subroutine fel_read_input (param_file, run, err_flag)

type (fel_run_struct), target :: run

! The namelist member names, as the input file writes them.
type (fel_global_struct) global
type (wavefront_init_struct) wavefront_init
type (fel_slicing_struct) slicing
type (fel_chamber_wake_init_struct) chamber_wake
type (fel_space_charge_input_struct) space_charge
type (beam_init_struct) beam_init
type (fel_resample_param_struct) resample
character(400) lat_file, beam_file, write_genesis_dist, write_openpmd_file
character(400) field_file(9)
logical quiet_start, shot_noise
character(16) load_mode
logical split_weights, swap_beam_xy, gen_test_weights, resample_split_weights
integer beamlet_size

character(*) param_file
logical err_flag

integer iu, ios
character(600) iomsg_text
character(*), parameter :: r_name = 'fel_read_input'

namelist / fel_params / lat_file, global, slicing, bmad_com, space_charge_com, chamber_wake, space_charge
namelist / fel_beam_init / beam_init, resample, beam_file, write_genesis_dist, &
                        write_openpmd_file, load_mode, quiet_start, beamlet_size, shot_noise, &
                        split_weights, swap_beam_xy, gen_test_weights, resample_split_weights
namelist / fel_wavefront_init / wavefront_init, field_file

!

err_flag = .true.

! Defaults: the declarations' own.

lat_file = run%lat_file
field_file = run%field_file
global = run%global
wavefront_init = run%winit
slicing = run%slicing
chamber_wake = run%chamber_wake
space_charge = run%space_charge
beam_init = run%beam_init
resample = run%resample
beam_file = run%bparam%beam_file
write_genesis_dist = run%bparam%write_genesis_dist
write_openpmd_file = run%bparam%write_openpmd_file
load_mode = run%bparam%load_mode
quiet_start = run%bparam%quiet_start
beamlet_size = run%bparam%beamlet_size
shot_noise = run%bparam%shot_noise
split_weights = run%bparam%split_weights
swap_beam_xy = run%bparam%swap_beam_xy
gen_test_weights = run%bparam%gen_test_weights
resample_split_weights = run%bparam%resample_split_weights

! The retired group is refused before anything is read.

if (has_retired_group(param_file)) return

open (newunit = iu, file = param_file, status = 'old', action = 'read', iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN INPUT FILE: ' // trim(param_file))
  return
endif

! Each group is optional, and an absent group keeps its defaults. A group that is
! present must parse.

if (group_present(iu, 'fel_params')) then
  rewind (iu)
  read (iu, nml = fel_params, iostat = ios, iomsg = iomsg_text)
  if (ios /= 0) then
    call out_io (s_error$, r_name, 'ERROR PARSING &fel_params IN: ' // trim(param_file), &
                                   trim(iomsg_text))
    close (iu)
    return
  endif
endif

if (group_present(iu, 'fel_beam_init')) then
  rewind (iu)
  read (iu, nml = fel_beam_init, iostat = ios, iomsg = iomsg_text)
  if (ios /= 0) then
    call out_io (s_error$, r_name, 'ERROR PARSING &fel_beam_init IN: ' // trim(param_file), &
                                   trim(iomsg_text))
    close (iu)
    return
  endif
endif

if (group_present(iu, 'fel_wavefront_init')) then
  rewind (iu)
  read (iu, nml = fel_wavefront_init, iostat = ios, iomsg = iomsg_text)
  if (ios /= 0) then
    call out_io (s_error$, r_name, 'ERROR PARSING &fel_wavefront_init IN: ' // trim(param_file), &
                                   trim(iomsg_text))
    close (iu)
    return
  endif
endif

close (iu)

run%lat_file = lat_file
run%field_file = field_file
run%global = global
run%winit = wavefront_init
run%slicing = slicing
run%chamber_wake = chamber_wake
run%space_charge = space_charge
run%beam_init = beam_init
run%resample = resample
run%bparam%beam_file = beam_file
run%bparam%write_genesis_dist = write_genesis_dist
run%bparam%write_openpmd_file = write_openpmd_file
run%bparam%load_mode = load_mode
run%bparam%quiet_start = quiet_start
run%bparam%beamlet_size = beamlet_size
run%bparam%shot_noise = shot_noise
run%bparam%split_weights = split_weights
run%bparam%swap_beam_xy = swap_beam_xy
run%bparam%gen_test_weights = gen_test_weights
run%bparam%resample_split_weights = resample_split_weights

err_flag = .false.

end subroutine fel_read_input

!------------------------------------------------------------------------------
!+
! Subroutine fel_write_resolved_input (run, iu)
!
! Routine to write the resolved inputs of run, every default made explicit, as the
! three namelist groups to an open unit. This is the provenance echo in the Meta/ group
! of the stats file, as Genesis's Meta group embeds the whole input file. An embedding
! program may also use it to save a configuration built in code.
!
! Input:
!   run -- fel_run_struct: The run state after parsing (resolved values).
!   iu  -- integer: Open unit to write the namelist echo to.
!
! Output:
!   None beyond the write: the three groups, every value resolved, in namelist syntax
!   (the Meta/input_echo content).
!-

subroutine fel_write_resolved_input (run, iu)

type (fel_run_struct), target :: run
integer iu

type (fel_global_struct) global
type (wavefront_init_struct) wavefront_init
type (fel_slicing_struct) slicing
type (fel_chamber_wake_init_struct) chamber_wake
type (fel_space_charge_input_struct) space_charge
type (beam_init_struct) beam_init
type (fel_resample_param_struct) resample
character(400) lat_file, beam_file, write_genesis_dist, write_openpmd_file
character(400) field_file(9)
logical quiet_start, shot_noise
character(16) load_mode
logical split_weights, swap_beam_xy, gen_test_weights, resample_split_weights
integer beamlet_size

namelist / fel_params / lat_file, global, slicing, bmad_com, space_charge_com, chamber_wake, space_charge
namelist / fel_beam_init / beam_init, resample, beam_file, write_genesis_dist, &
                        write_openpmd_file, load_mode, quiet_start, beamlet_size, shot_noise, &
                        split_weights, swap_beam_xy, gen_test_weights, resample_split_weights
namelist / fel_wavefront_init / wavefront_init, field_file

!

lat_file = run%lat_file
field_file = run%field_file
global = run%global
wavefront_init = run%winit
slicing = run%slicing
chamber_wake = run%chamber_wake
space_charge = run%space_charge
beam_init = run%beam_init
resample = run%resample
beam_file = run%bparam%beam_file
write_genesis_dist = run%bparam%write_genesis_dist
write_openpmd_file = run%bparam%write_openpmd_file
load_mode = run%bparam%load_mode
quiet_start = run%bparam%quiet_start
beamlet_size = run%bparam%beamlet_size
shot_noise = run%bparam%shot_noise
split_weights = run%bparam%split_weights
swap_beam_xy = run%bparam%swap_beam_xy
gen_test_weights = run%bparam%gen_test_weights
resample_split_weights = run%bparam%resample_split_weights

write (iu, nml = fel_params)
write (iu, nml = fel_beam_init)
write (iu, nml = fel_wavefront_init)

end subroutine fel_write_resolved_input

!------------------------------------------------------------------------------
!+
! Function group_present (iu, group_name) result (found)
!
! Routine to find whether the namelist group &group_name is in the open file. The scan is
! textual, since the standard gives no portable iostat for an absent group. Leading
! blanks and case are ignored.
!
! Input:
!   iu         -- integer: Open unit positioned at the file start.
!   group_name -- character(*): Namelist group name to look for, without the &.
!
! Output:
!   found      -- logical: True when the group appears in the file.
!-

function group_present (iu, group_name) result (found)

integer iu, ios
character(*) group_name
character(400) line
logical found

!

found = .false.
rewind (iu)
do
  read (iu, '(a)', iostat = ios) line
  if (ios /= 0) exit
  line = adjustl(line)
  if (line(1:1) /= '&') cycle
  call downcase_string (line)
  if (line(2:len_trim(group_name)+1) == trim(group_name) .and. &
      (len_trim(line) == len_trim(group_name)+1 .or. &
       line(len_trim(group_name)+2:len_trim(group_name)+2) == ' ')) then
    found = .true.
    exit
  endif
enddo
rewind (iu)

end function group_present

!------------------------------------------------------------------------------
!+
! Function has_retired_group (param_file) result (found)
!
! Routine to refuse retired input names, each with the name that replaced it. So
! migration is a mechanical edit of the input file. Refused are the flat
! &fel_track_params group, the window parameters that moved from wavefront_init to
! slicing, dist_file and use_beam_init.
!
! Input:
!   param_file -- character(*): Input file name.
!
! Output:
!   found      -- logical: True when a retired name is present. The caller then refuses
!                   the run, and the message names what replaced it.
!-

function has_retired_group (param_file) result (found)

character(*) param_file
logical found

integer iu, ios
character(400) line
character(*), parameter :: r_name = 'fel_read_input'

! The scan is textual. Fortran ignores an unknown namelist group in silence. A retired
! component of a live group is a parse error, and its message names the component but
! not what to do. A deck written against either would otherwise run on defaults or stop
! without advice.

found = .false.
open (newunit = iu, file = param_file, status = 'old', action = 'read', iostat = ios)
if (ios /= 0) return

do
  read (iu, '(a)', iostat = ios) line
  if (ios /= 0) exit
  line = adjustl(line)
  call downcase_string (line)
  if (line(1:17) == '&fel_track_params') then
    found = .true.
    call out_io (s_error$, r_name, '&fel_track_params IS NOT AN INPUT GROUP.', &
      'THE THREE GROUPS ARE &fel_params, &fel_beam_init AND &fel_wavefront_init.', &
      'SEE docs/input-reference.md.')
    exit
  endif

  ! The window belongs to the FEL interaction and not to the field. So it moved from
  ! wavefront_init to slicing when the wavelengths per slice became an integer.

  if (index(line, 'wavefront_init%window_sample') /= 0) then
    found = .true.
    call out_io (s_error$, r_name, 'wavefront_init%window_sample IS RETIRED.', &
      'POSSIBLE SOLUTION: SET slicing%n_wavelength IN &fel_params.')
    exit
  endif
  if (index(line, 'wavefront_init%window_length') /= 0) then
    found = .true.
    call out_io (s_error$, r_name, 'wavefront_init%window_length IS RETIRED.', &
      'POSSIBLE SOLUTION: SET slicing%window_length IN &fel_params.')
    exit
  endif

  ! A bunch to slice belongs to beam_init, read or generated, and Bmad's position_file
  ! reads it. The flag that resampled a generated bunch is now resample%use_beam_init.

  if (index(line, 'dist_file') /= 0 .and. index(line, 'resample%') == 0) then
    found = .true.
    call out_io (s_error$, r_name, 'dist_file IS RETIRED.', &
      'POSSIBLE SOLUTION: SET beam_init%position_file. THE BUNCH IS THEN SLICED BY load_mode.')
    exit
  endif
  if (index(line, 'use_beam_init') /= 0 .and. index(line, 'resample%') == 0) then
    found = .true.
    call out_io (s_error$, r_name, 'use_beam_init IS RETIRED.', &
      'A GENERATED BUNCH TAKES THE ANALYTIC LOADER. TO RESAMPLE ONE FOR A CHECK, SET resample%use_beam_init.')
    exit
  endif
enddo

close (iu)

end function has_retired_group

end module fel_input_mod
