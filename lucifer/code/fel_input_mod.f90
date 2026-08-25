!+
! Module fel_input_mod
!
! The namelist layer of the FEL tracker, quarantined the way Tao quarantines its
! tao_init_* files: three groups, all read from ONE input file, each filling the
! structs of fel_struct directly (Tao's &tao_params pattern: set global%out_root,
! bmad_com%radiation_damping_on, wake%radius, sc%nz, beam_init%n_particle,
! wavefront_init%lambda0 by component). This module is itself library: an embedding
! program may reuse the parsing, or skip it and fill the structs in code.
!
!   &fel_params          lat_file, global, bmad_com, space_charge_com, wake, sc,
!                        write_wake_kernels
!   &fel_beam_init       beam_init, imp, beam_file, dist_file, write_dist_file,
!                        write_opmd_file, use_beam_init, nbins, shotnoise,
!                        split_weights, swap_beam_xy, gen_test_weights,
!                        imp_split_weights
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
type (fel_wake_init_struct) wake
type (fel_efield_struct) sc
type (beam_init_struct) beam_init
type (fel_import_param_struct) imp
character(400) lat_file, beam_file, dist_file, write_dist_file, write_opmd_file
character(400) write_wake_kernels
character(400) field_file(9)
logical use_beam_init, shotnoise
logical split_weights, swap_beam_xy, gen_test_weights, imp_split_weights
integer nbins

character(*) param_file
logical err_flag

integer iu, ios
character(600) iomsg_text
character(*), parameter :: r_name = 'fel_read_input'

namelist / fel_params / lat_file, global, bmad_com, space_charge_com, wake, sc, &
                        write_wake_kernels
namelist / fel_beam_init / beam_init, imp, beam_file, dist_file, write_dist_file, &
                        write_opmd_file, use_beam_init, nbins, shotnoise, &
                        split_weights, swap_beam_xy, gen_test_weights, imp_split_weights
namelist / fel_wavefront_init / wavefront_init, field_file

!

err_flag = .true.

! Defaults: the declarations' own.

lat_file = run%lat_file
field_file = run%field_file
global = run%global
wavefront_init = run%winit
wake = run%wake_init
sc = run%sc_init
beam_init = run%beam_init
imp = run%imp
beam_file = run%bparam%beam_file
dist_file = run%bparam%dist_file
write_dist_file = run%bparam%write_dist_file
write_opmd_file = run%bparam%write_opmd_file
use_beam_init = run%bparam%use_beam_init
nbins = run%bparam%nbins
shotnoise = run%bparam%shotnoise
split_weights = run%bparam%split_weights
swap_beam_xy = run%bparam%swap_beam_xy
gen_test_weights = run%bparam%gen_test_weights
imp_split_weights = run%bparam%imp_split_weights
write_wake_kernels = run%wake_init%write_kernels

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
run%wake_init = wake
run%wake_init%write_kernels = write_wake_kernels
run%sc_init = sc
run%beam_init = beam_init
run%imp = imp
run%bparam%beam_file = beam_file
run%bparam%dist_file = dist_file
run%bparam%write_dist_file = write_dist_file
run%bparam%write_opmd_file = write_opmd_file
run%bparam%use_beam_init = use_beam_init
run%bparam%nbins = nbins
run%bparam%shotnoise = shotnoise
run%bparam%split_weights = split_weights
run%bparam%swap_beam_xy = swap_beam_xy
run%bparam%gen_test_weights = gen_test_weights
run%bparam%imp_split_weights = imp_split_weights

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
type (fel_wake_init_struct) wake
type (fel_efield_struct) sc
type (beam_init_struct) beam_init
type (fel_import_param_struct) imp
character(400) lat_file, beam_file, dist_file, write_dist_file, write_opmd_file
character(400) write_wake_kernels
character(400) field_file(9)
logical use_beam_init, shotnoise
logical split_weights, swap_beam_xy, gen_test_weights, imp_split_weights
integer nbins

namelist / fel_params / lat_file, global, bmad_com, space_charge_com, wake, sc, &
                        write_wake_kernels
namelist / fel_beam_init / beam_init, imp, beam_file, dist_file, write_dist_file, &
                        write_opmd_file, use_beam_init, nbins, shotnoise, &
                        split_weights, swap_beam_xy, gen_test_weights, imp_split_weights
namelist / fel_wavefront_init / wavefront_init, field_file

!

lat_file = run%lat_file
field_file = run%field_file
global = run%global
wavefront_init = run%winit
wake = run%wake_init
sc = run%sc_init
beam_init = run%beam_init
imp = run%imp
beam_file = run%bparam%beam_file
dist_file = run%bparam%dist_file
write_dist_file = run%bparam%write_dist_file
write_opmd_file = run%bparam%write_opmd_file
use_beam_init = run%bparam%use_beam_init
nbins = run%bparam%nbins
shotnoise = run%bparam%shotnoise
split_weights = run%bparam%split_weights
swap_beam_xy = run%bparam%swap_beam_xy
gen_test_weights = run%bparam%gen_test_weights
imp_split_weights = run%bparam%imp_split_weights
write_wake_kernels = run%wake_init%write_kernels

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

integer iu, ios, i, j
character(400) line
character(60) key
logical in_group
character(*), parameter :: r_name = 'fel_read_input'

! old-name -> new-home table. Names not listed moved verbatim under the shown group.

character(56), parameter :: moved(2, 44) = reshape([ character(56) :: &
  'lat_file',           '&fel_params lat_file',                        &
  'out_root',           '&fel_params global%out_root',                 &
  'interlude_model',    '&fel_params global%interlude_model',          &
  'wavefront_format',   '&fel_params global%wavefront_formats (list)', &
  'write_diag',         '&fel_params global%write_diag',               &
  'write_initial',      '&fel_params global%write_initial',            &
  'load_only',          '&fel_params global%load_only',                &
  'keep_escaped_field', '&fel_params global%keep_escaped_field',       &
  'dump_beam_at',       '&fel_params global%dump_beam_at',             &
  'dump_field_at',      '&fel_params global%dump_field_at',            &
  'ran_seed',           '&fel_params global%ran_seed',                 &
  'migrate',            '&fel_params global%migrate',                  &
  'migrate_check',      '&fel_params global%migrate_check',            &
  'reference_run',      '&fel_params global%reference_run',            &
  'radiation_damping',  '&fel_params bmad_com%radiation_damping_on',   &
  'radiation_fluctuations', '&fel_params bmad_com%radiation_fluctuations_on', &
  'wake_on',            '&fel_params wake%on',                         &
  'wake_loss',          '&fel_params wake%loss',                       &
  'wake_radius',        '&fel_params wake%radius',                     &
  'wake_conductivity',  '&fel_params wake%conductivity',               &
  'wake_relaxation',    '&fel_params wake%relaxation',                 &
  'wake_roundpipe',     '&fel_params wake%roundpipe',                  &
  'wake_material',      '&fel_params wake%material',                   &
  'wake_gap',           '&fel_params wake%gap',                        &
  'wake_lgap',          '&fel_params wake%lgap',                       &
  'wake_hrough',        '&fel_params wake%hrough',                     &
  'wake_lrough',        '&fel_params wake%lrough',                     &
  'sc_rmax',            '&fel_params sc%rmax',                         &
  'sc_ngrid',           '&fel_params sc%ngrid',                        &
  'sc_nz',              '&fel_params sc%nz',                           &
  'sc_nphi',            '&fel_params sc%nphi',                         &
  'sc_longrange',       '&fel_params sc%longrange',                    &
  'write_wake_kernels', '&fel_params write_wake_kernels',              &
  'lambda0',            '&fel_wavefront_init wavefront_init%lambda0',  &
  'window_length',      '&fel_wavefront_init wavefront_init%window_length', &
  'window_sample',      '&fel_wavefront_init wavefront_init%window_sample', &
  'grid_n_pts',         '&fel_wavefront_init wavefront_init%grid_n_pts', &
  'grid_half_width',    '&fel_wavefront_init wavefront_init%grid_half_width', &
  'seed_power',         '&fel_wavefront_init wavefront_init%seed_power', &
  'seed_waist_size',    '&fel_wavefront_init wavefront_init%seed_waist_size', &
  'seed_polarization',  '&fel_wavefront_init wavefront_init%seed_polarization', &
  'harmonics',          '&fel_wavefront_init wavefront_init%harmonics', &
  'field_file',         '&fel_wavefront_init field_file',              &
  'beam_file',          '&fel_beam_init beam_file'], [2, 44])

!

found = .false.
open (newunit = iu, file = param_file, status = 'old', action = 'read', iostat = ios)
if (ios /= 0) return

in_group = .false.
do
  read (iu, '(a)', iostat = ios) line
  if (ios /= 0) exit
  line = adjustl(line)
  call downcase_string (line)
  if (line(1:17) == '&fel_track_params') then
    in_group = .true.
    found = .true.
    call out_io (s_error$, r_name, '&fel_track_params IS RETIRED. Its parameters moved to the', &
      'three groups &fel_params, &fel_beam_init and &fel_wavefront_init. In this file:')
    cycle
  endif
  if (.not. in_group) cycle
  if (line(1:1) == '/' .or. line(1:4) == '&end') exit
  if (line(1:1) == '!' .or. line == '') cycle
  i = index(line, '=')
  if (i <= 1) cycle
  key = adjustl(line(1:i-1))
  i = index(key, '(')                     ! field_file(2) -> field_file
  if (i > 1) key = key(1:i-1)
  i = index(key, '%')                     ! beam_init%sig_z -> beam_init
  if (i > 1) key = key(1:i-1)
  do j = 1, size(moved, 2)
    if (key == moved(1, j)) then
      call out_io (s_blank$, r_name, '  ' // trim(key) // ' -> ' // trim(moved(2, j)))
      exit
    endif
  enddo
  if (key == 'beam_init' .or. key == 'imp' .or. key == 'dist_file' .or. &
      key == 'write_dist_file' .or. key == 'write_opmd_file' .or. key == 'use_beam_init' .or. &
      key == 'nbins' .or. key == 'shotnoise' .or. key == 'split_weights' .or. &
      key == 'swap_beam_xy' .or. key == 'gen_test_weights' .or. key == 'imp_split_weights') then
    call out_io (s_blank$, r_name, '  ' // trim(key) // ' -> &fel_beam_init ' // trim(key))
  endif
enddo
close (iu)

end function has_retired_group

end module fel_input_mod
