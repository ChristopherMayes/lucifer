!+
! Module fel_record_mod
!
! The numbers the printout shows, kept as values for the stats file.
!
! The header and the footer print a number and hand the same number here, under a group
! and a name, with its unit and a description. fel_record_write puts every one of them in
! the stats file as a scalar dataset, setup/ for what setup worked out before tracking,
! run/ for what the run produced and timing/ for the wall clock, which is the one group
! that differs between two runs of one input, and the files the run wrote as the attributes
! files and file_bytes of run/. A program that wants a number the printout shows reads it there, so the
! printout's layout is free to change.
!
! The entries are held in fel_record_com between the call that hands one in and the
! write. fel_record_write may be called more than once on one file, and each call writes
! the entries handed in since the last.
!-

module fel_record_mod

use fel_h5_mod

implicit none

private
public :: fel_record_real, fel_record_int, fel_record_flag, fel_record_str, fel_record_file
public :: fel_record_write, fel_record_com, fel_record_common_struct, fel_record_entry_struct

integer, parameter :: fel_record_real$ = 1, fel_record_int$ = 2, fel_record_flag$ = 3, fel_record_str$ = 4

type fel_record_entry_struct
  character(16) :: group = ''            ! setup, run or timing.
  character(48) :: name = ''             ! The dataset's name.
  character(16) :: unit = ''             ! Its unit, 1 for a pure number.
  character(300) :: descrip = ''         ! Its description.
  integer :: kind = 0                    ! fel_record_real$, _int$, _flag$ or _str$.
  real(rp) :: r = 0                      ! A real value.
  integer :: i = 0                       ! An integer value.
  logical :: flag = .false.              ! A logical value.
  character(300) :: s = ''               ! A string value.
  logical :: written = .false.           ! Already in the file.
end type

type fel_record_common_struct
  type (fel_record_entry_struct), allocatable :: entry(:)   ! The entries, 1 to n.
  integer :: n = 0                                          ! Entries held.
  character(300), allocatable :: file(:)                    ! Files the run wrote, 1 to n_file.
  real(rp), allocatable :: file_bytes(:)                    ! Their sizes (bytes).
  integer :: n_file = 0                                     ! Files held.
end type

type (fel_record_common_struct), save, target :: fel_record_com

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function new_entry (group, name, unit, descrip, kind) result (ix)
!
! Routine to make room for one entry and fill what every kind shares. An entry of the
! same group and name is replaced, so a number printed twice is held once.
!
! Input:
!   group   -- character(*): setup, run or timing.
!   name    -- character(*): The dataset's name.
!   unit    -- character(*): Its unit.
!   descrip -- character(*): Its description.
!   kind    -- integer: The kind of value.
!
! Output:
!   ix      -- integer: The entry's index in fel_record_com%entry.
!-

function new_entry (group, name, unit, descrip, kind) result (ix)

type (fel_record_entry_struct), allocatable :: more(:)
character(*) group, name, unit, descrip
integer kind, ix

!

if (.not. allocated(fel_record_com%entry)) allocate (fel_record_com%entry(64))

do ix = 1, fel_record_com%n
  if (fel_record_com%entry(ix)%group == group .and. fel_record_com%entry(ix)%name == name .and. &
        .not. fel_record_com%entry(ix)%written) exit
enddo

if (ix > fel_record_com%n) then
  if (fel_record_com%n == size(fel_record_com%entry)) then
    allocate (more(2 * fel_record_com%n))
    more(1:fel_record_com%n) = fel_record_com%entry
    call move_alloc (more, fel_record_com%entry)
  endif
  fel_record_com%n = fel_record_com%n + 1
  ix = fel_record_com%n
endif

fel_record_com%entry(ix) = fel_record_entry_struct(group = group, name = name, unit = unit, &
                                                    descrip = descrip, kind = kind)

end function new_entry

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_record_real (group, name, unit, descrip, value)
!
! Routine to hold one real number for the stats file.
!
! Input:
!   group   -- character(*): setup, run or timing.
!   name    -- character(*): The dataset's name.
!   unit    -- character(*): Its unit, 1 for a pure number.
!   descrip -- character(*): Its description.
!   value   -- real(rp): The number, in the unit named.
!-

subroutine fel_record_real (group, name, unit, descrip, value)

character(*) group, name, unit, descrip
real(rp) value
integer ix

!

ix = new_entry(group, name, unit, descrip, fel_record_real$)
fel_record_com%entry(ix)%r = value

end subroutine fel_record_real

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_record_int (group, name, descrip, value)
!
! Routine to hold one count for the stats file.
!
! Input:
!   group   -- character(*): setup, run or timing.
!   name    -- character(*): The dataset's name.
!   descrip -- character(*): Its description.
!   value   -- integer: The count.
!-

subroutine fel_record_int (group, name, descrip, value)

character(*) group, name, descrip
integer value, ix

!

ix = new_entry(group, name, '1', descrip, fel_record_int$)
fel_record_com%entry(ix)%i = value

end subroutine fel_record_int

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_record_flag (group, name, descrip, value)
!
! Routine to hold one logical for the stats file.
!
! Input:
!   group   -- character(*): setup, run or timing.
!   name    -- character(*): The dataset's name.
!   descrip -- character(*): Its description.
!   value   -- logical: The value.
!-

subroutine fel_record_flag (group, name, descrip, value)

character(*) group, name, descrip
logical value
integer ix

!

ix = new_entry(group, name, '1', descrip, fel_record_flag$)
fel_record_com%entry(ix)%flag = value

end subroutine fel_record_flag

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_record_str (group, name, descrip, value)
!
! Routine to hold one string for the stats file.
!
! Input:
!   group   -- character(*): setup, run or timing.
!   name    -- character(*): The dataset's name.
!   descrip -- character(*): Its description.
!   value   -- character(*): The string.
!-

subroutine fel_record_str (group, name, descrip, value)

character(*) group, name, descrip, value
integer ix

!

ix = new_entry(group, name, '', descrip, fel_record_str$)
fel_record_com%entry(ix)%s = value

end subroutine fel_record_str

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_record_file (name, bytes)
!
! Routine to hold the name of one file the run wrote, for the attribute files of run/. A series of frames
! is one name in shell brace syntax, as the printout lists it.
!
! Input:
!   name  -- character(*): The file's name, relative to the run's directory.
!   bytes -- real(rp): Its size (bytes), the sum over a series.
!-

subroutine fel_record_file (name, bytes)

character(*) name
real(rp) bytes
character(300), allocatable :: more(:)
real(rp), allocatable :: more_bytes(:)

!

if (.not. allocated(fel_record_com%file)) allocate (fel_record_com%file(32), fel_record_com%file_bytes(32))
if (fel_record_com%n_file == size(fel_record_com%file)) then
  allocate (more(2 * fel_record_com%n_file), more_bytes(2 * fel_record_com%n_file))
  more(1:fel_record_com%n_file) = fel_record_com%file
  more_bytes(1:fel_record_com%n_file) = fel_record_com%file_bytes
  call move_alloc (more, fel_record_com%file)
  call move_alloc (more_bytes, fel_record_com%file_bytes)
endif
fel_record_com%n_file = fel_record_com%n_file + 1
fel_record_com%file(fel_record_com%n_file) = name
fel_record_com%file_bytes(fel_record_com%n_file) = bytes

end subroutine fel_record_file

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_record_write (file_name, with_files, err_flag)
!
! Routine to write the entries held, and with with_files the list of files, into a stats
! file that is already there. The groups setup/, run/ and timing/ are made where the
! file has none. An entry is written once: a second call writes what was handed in since.
!
! Input:
!   file_name  -- character(*): The stats file.
!   with_files -- logical: Write the attributes files and file_bytes of run/ too. The
!                   stats file's own size is its size before this call.
!
! Output:
!   err_flag   -- logical: Set true if the file could not be opened or written.
!-

subroutine fel_record_write (file_name, with_files, err_flag)

type (fel_record_entry_struct), pointer :: en
character(*) file_name
integer(hid_t) f_id, g_id
integer ig, ix, h5e
logical with_files, err_flag, err, there, ferr
character(16), parameter :: groups(3) = [character(16):: 'setup', 'run', 'timing']
character(90), parameter :: group_descrip(3) = [character(90):: &
      'What setup worked out from the input before tracking, one scalar each.', &
      'What the run produced, one scalar each.', &
      'The wall clock of the run by phase, its threads, and the estimate made before tracking.']
character(*), parameter :: r_name = 'fel_record_write'

!

err_flag = .true.
call hdf5_open_file (file_name, 'APPEND', f_id, err)
if (err) then
  call out_io (s_error$, r_name, 'Could not reopen the stats file ' // trim(file_name) // '.')
  return
endif

err = .false.
do ig = 1, size(groups)
  if (.not. (any(fel_record_com%entry(1:fel_record_com%n)%group == groups(ig) .and. &
             .not. fel_record_com%entry(1:fel_record_com%n)%written) .or. &
             (with_files .and. groups(ig) == 'run'))) cycle

  call H5Lexists_f (f_id, trim(groups(ig)), there, h5e)
  if (there) then
    call H5Gopen_f (f_id, trim(groups(ig)), g_id, h5e)
  else
    call H5Gcreate_f (f_id, trim(groups(ig)), g_id, h5e)
    call H5LTset_attribute_string_f (f_id, trim(groups(ig)), 'kind', 'run', h5e)
    call H5LTset_attribute_string_f (f_id, trim(groups(ig)), 'description', trim(group_descrip(ig)), h5e)
  endif
  err = err .or. (h5e < 0)

  do ix = 1, fel_record_com%n
    en => fel_record_com%entry(ix)
    if (en%group /= groups(ig) .or. en%written) cycle
    select case (en%kind)
    case (fel_record_real$)
      call fel_h5_real (g_id, trim(en%name), trim(en%unit), label_of(en%name), trim(en%descrip), '', en%r, err)
    case (fel_record_int$)
      call fel_h5_int (g_id, trim(en%name), '1', label_of(en%name), trim(en%descrip), '', en%i, err)
    case (fel_record_flag$)
      call fel_h5_flag (g_id, trim(en%name), label_of(en%name), trim(en%descrip), en%flag, err)
    case (fel_record_str$)
      call fel_h5_text (g_id, trim(en%name), label_of(en%name), trim(en%descrip), trim(en%s), err)
    end select
    en%written = .true.
  enddo

  if (with_files .and. groups(ig) == 'run' .and. fel_record_com%n_file > 0) then
    ! Attributes of run/, as the stats specification has a list with no axis.
    call hdf5_write_attribute_string_rank1 (g_id, 'files', fel_record_com%file(1:fel_record_com%n_file), ferr)
    err = err .or. ferr
    call hdf5_write_attribute_real_rank1 (g_id, 'file_bytes', fel_record_com%file_bytes(1:fel_record_com%n_file), ferr)
    err = err .or. ferr
  endif

  call H5Gclose_f (g_id, h5e)
enddo

call H5Fclose_f (f_id, h5e)
err_flag = err

!-------------------------------------------------------
contains

!+
! Function label_of (name) result (label)
!
! Routine to make a dataset's label from its name, the underscores as spaces.
!-

function label_of (name) result (label)

character(*) name
character(:), allocatable :: label
integer i

!

label = trim(name)
do i = 1, len(label)
  if (label(i:i) == '_') label(i:i) = ' '
enddo

end function label_of

end subroutine fel_record_write

end module fel_record_mod
