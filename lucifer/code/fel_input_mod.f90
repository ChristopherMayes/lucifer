!+
! Module fel_input_mod
!
! The namelist layer of the FEL tracker, quarantined the way Tao quarantines its
! tao_init_* files: three groups, all read from ONE input file, each filling the
! structs of fel_struct directly (Tao's &tao_params pattern: set global%out_root,
! bmad_com%radiation_damping_on, chamber_wake%radius, space_charge%nz, beam_init%n_particle,
! wavefront_init%lambda0 by component). This module is itself library: an embedding
! program may reuse the parsing, or skip it and fill the structs in code.
!
!   &fel_params          lat_file, global, bmad_com, space_charge_com, wake, sc,
!                        chamber_wake%write_kernels
!   &fel_beam_init       beam_init, imp, beam_file, dist_file, write_genesis_dist,
!                        write_openpmd_file, use_beam_init, beamlet_size, shot_noise,
!                        split_weights, swap_beam_xy, gen_test_weights,
!                        resample_split_weights
!   &fel_wavefront_init  wavefront_init, field_file
!
! A group that is absent keeps its defaults. A group that is present must parse.
! The retired flat &fel_track_params group is refused BY NAME: the error lists each
! parameter found in it and the group it moved to.
!-

module fel_input_mod

use fel_struct

implicit none

contains

!------------------------------------------------------------------------------
!+
! Subroutine fel_read_input (param_file, run, err_flag)
!
! Read the three namelist groups of param_file into run's input structs. bmad_com and
! space_charge_com are set directly by the namelist (they are Bmad's own globals).
!
! Input:
!   param_file   -- character(*): the input file.
!
! Output:
!   run                -- fel_run_struct: input structs filled (global, winit,
!                           wake_init, sc_init, beam_init, imp, bparam, lat_file,
!                           field_file). Nothing else is touched.
!   err_flag           -- logical: set on any parse error or on the retired group.
!-

subroutine fel_read_input (param_file, run, err_flag)

type (fel_run_struct), target :: run

! The namelist member names, as the input file writes them.
type (fel_global_struct) global
type (wavefront_init_struct) wavefront_init
type (fel_chamber_wake_init_struct) chamber_wake
type (fel_space_charge_struct) space_charge
type (beam_init_struct) beam_init
type (fel_resample_param_struct) resample
character(400) lat_file, beam_file, dist_file, write_genesis_dist, write_openpmd_file
character(400) field_file(9)
logical use_beam_init, shot_noise
logical split_weights, swap_beam_xy, gen_test_weights, resample_split_weights
integer beamlet_size

character(*) param_file
logical err_flag

integer iu, ios
character(600) iomsg_text
character(*), parameter :: r_name = 'fel_read_input'

namelist / fel_params / lat_file, global, bmad_com, space_charge_com, chamber_wake, space_charge
namelist / fel_beam_init / beam_init, resample, beam_file, dist_file, write_genesis_dist, &
                        write_openpmd_file, use_beam_init, beamlet_size, shot_noise, &
                        split_weights, swap_beam_xy, gen_test_weights, resample_split_weights
namelist / fel_wavefront_init / wavefront_init, field_file

!

err_flag = .true.

! Defaults: the declarations' own.

lat_file = run%lat_file
field_file = run%field_file
global = run%global
wavefront_init = run%winit
chamber_wake = run%chamber_wake
space_charge = run%space_charge
beam_init = run%beam_init
resample = run%resample
beam_file = run%bparam%beam_file
dist_file = run%bparam%dist_file
write_genesis_dist = run%bparam%write_genesis_dist
write_openpmd_file = run%bparam%write_openpmd_file
use_beam_init = run%bparam%use_beam_init
beamlet_size = run%bparam%beamlet_size
shot_noise = run%bparam%shot_noise
split_weights = run%bparam%split_weights
swap_beam_xy = run%bparam%swap_beam_xy
gen_test_weights = run%bparam%gen_test_weights
resample_split_weights = run%bparam%resample_split_weights

! The retired group is refused by name before anything is read.

if (has_retired_group(param_file)) return

open (newunit = iu, file = param_file, status = 'old', action = 'read', iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN INPUT FILE: ' // trim(param_file))
  return
endif

! Each group is optional (absent = defaults) but must parse when present.

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
run%chamber_wake = chamber_wake
run%space_charge = space_charge
run%beam_init = beam_init
run%resample = resample
run%bparam%beam_file = beam_file
run%bparam%dist_file = dist_file
run%bparam%write_genesis_dist = write_genesis_dist
run%bparam%write_openpmd_file = write_openpmd_file
run%bparam%use_beam_init = use_beam_init
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
! Write run's RESOLVED inputs (every default made explicit) as the three
! namelist groups, to an open unit. The stats file's Meta/ provenance echo
! (Genesis parity: its Meta group embeds the entire input file). An embedding
! program may also use it to persist a configuration built in code.
!
! Input:
!   run -- fel_run_struct: The run state after parsing (resolved values).
!   iu  -- integer: Open unit to write the namelist echo to.
!
! Output:
!   None beyond the write: the three groups, every value resolved, in namelist
!   syntax (the Meta/input_echo content).
!-

subroutine fel_write_resolved_input (run, iu)

type (fel_run_struct), target :: run
integer iu

type (fel_global_struct) global
type (wavefront_init_struct) wavefront_init
type (fel_chamber_wake_init_struct) chamber_wake
type (fel_space_charge_struct) space_charge
type (beam_init_struct) beam_init
type (fel_resample_param_struct) resample
character(400) lat_file, beam_file, dist_file, write_genesis_dist, write_openpmd_file
character(400) field_file(9)
logical use_beam_init, shot_noise
logical split_weights, swap_beam_xy, gen_test_weights, resample_split_weights
integer beamlet_size

namelist / fel_params / lat_file, global, bmad_com, space_charge_com, chamber_wake, space_charge
namelist / fel_beam_init / beam_init, resample, beam_file, dist_file, write_genesis_dist, &
                        write_openpmd_file, use_beam_init, beamlet_size, shot_noise, &
                        split_weights, swap_beam_xy, gen_test_weights, resample_split_weights
namelist / fel_wavefront_init / wavefront_init, field_file

!

lat_file = run%lat_file
field_file = run%field_file
global = run%global
wavefront_init = run%winit
chamber_wake = run%chamber_wake
space_charge = run%space_charge
beam_init = run%beam_init
resample = run%resample
beam_file = run%bparam%beam_file
dist_file = run%bparam%dist_file
write_genesis_dist = run%bparam%write_genesis_dist
write_openpmd_file = run%bparam%write_openpmd_file
use_beam_init = run%bparam%use_beam_init
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
! Is the namelist group &group_name present in the (open) file? A textual scan
! (the standard gives no portable "group absent" iostat), tolerant of leading blanks
! and case.
!
! Input:
!   iu         -- integer: Open unit positioned at the file start.
!   group_name -- character(*): Namelist group name to look for (with the &).
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
! Refuse the retired flat &fel_track_params group BY NAME: list every parameter set
! in it together with the group and name it moved to, so migration is a mechanical
! edit of the input file.
!
! Input:
!   param_file -- character(*): Input file name.
!
! Output:
!   found      -- logical: True when the retired &fel_track_params group is present
!                   (the caller then refuses with the parameter mapping table).
!-

function has_retired_group (param_file) result (found)

character(*) param_file
logical found

integer iu, ios
character(400) line
character(*), parameter :: r_name = 'fel_read_input'

! A textual scan, because Fortran ignores an unknown namelist group in silence: a deck
! written against the retired flat group would otherwise run on defaults.

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
      'SEE doc/input-reference.md.')
    exit
  endif
enddo

close (iu)

end function has_retired_group

end module fel_input_mod
