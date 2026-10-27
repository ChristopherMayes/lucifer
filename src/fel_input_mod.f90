!+
! Module fel_input_mod
!
! The input reader of the FEL tracker, kept apart the way Tao keeps its tao_init_*
! files. An input file is one TOML file, read through the schema (fel_schema_mod) into the
! input structures of fel_struct. This module is library code. An embedding program
! may use the reader, or skip it and fill the structures in code.
!
! Every key the input file can set is in the schema's table, and a key the table does not
! name is refused with the keys of its table, so a misspelling cannot pass as a
! default. A key in the wrong table is refused with the table it belongs in. The lattice
! is required, and everything else keeps the declarations' defaults when absent.
!
! A key can also be set on the command line, after the input file, as table.key=value with the
! value in TOML syntax, lucifer input.toml beam.n_particle=4096 run.out_root=s7. A bare
! value for a string key is quoted, since the shell strips the quotes. The overrides are
! one TOML text parsed in memory and read through the same walk after the input file, so an
! override is refused as an input line is.
!-

module fel_input_mod

use fel_struct
use fel_schema_mod
use tomlf, only: toml_table, toml_key, toml_error, toml_load, toml_loads, get_value

implicit none

contains

!------------------------------------------------------------------------------
!+
! Subroutine fel_read_input (input_file, run, err_flag, overrides)
!
! Routine to read a TOML input file into the input structures of run, and then the overrides
! given on the command line. The Bmad switches of the [bmad] table set bmad_com and
! space_charge_com directly, since they are Bmad's own globals.
!
! Input:
!   input_file  -- character(*): The input file.
!   overrides  -- character(*), optional: The arguments after the input file, each
!                   table.key=value or lattice=value with the value in TOML syntax.
!
! Output:
!   run        -- fel_run_struct: Input structures filled (global, winit, slicing,
!                   chamber_wake, space_charge, beam_init, resample, bparam, lat_file,
!                   field_file). Nothing else is touched.
!   err_flag   -- logical: Set when the file does not parse, a key is unknown or in the
!                   wrong table, a value has the wrong type or is outside its bounds, the
!                   lattice is missing, or an override is not table.key=value.
!-

subroutine fel_read_input (input_file, run, err_flag, overrides)

type (fel_run_struct), target :: run
type (fel_key_struct), allocatable :: keys(:)
type (toml_table), allocatable :: root
type (toml_table), pointer :: tab, sub
type (toml_key), allocatable :: tkeys(:), ckeys(:)
type (toml_error), allocatable :: terr
character(*) input_file
character(*), optional :: overrides(:)
character(:), allocatable :: text, origin
logical err_flag, found_lattice, bad
integer i, j
character(*), parameter :: r_name = 'fel_read_input'

!

err_flag = .true.
found_lattice = .false.
call fel_schema_table (keys)

call toml_load (root, trim(input_file), error = terr)
if (allocated(terr)) then
  call out_io (s_error$, r_name, 'CANNOT READ THE INPUT FILE ' // trim(input_file) // ':', trim(terr%message))
  return
endif

origin = 'THE INPUT FILE'
if (.not. read_root()) return

if (.not. found_lattice) then
  call out_io (s_error$, r_name, 'THE INPUT FILE ' // trim(input_file) // ' NAMES NO LATTICE. ' // &
        'The top-level key lattice is required.')
  return
endif

! The overrides, as one TOML text of dotted keys, beam.n_particle = 4096, which TOML
! reads into the table beam. A heading per override, [beam], would be refused on the
! second key of one table, since a table can be headed once.

if (present(overrides)) then
  if (size(overrides) > 0) then
    text = ''
    do i = 1, size(overrides)
      if (.not. override_text(overrides(i))) return
    enddo
    deallocate (root)
    call toml_loads (root, text, error = terr)
    if (allocated(terr)) then
      call out_io (s_error$, r_name, 'CANNOT READ THE OVERRIDES ON THE COMMAND LINE:', trim(terr%message))
      return
    endif
    origin = 'THE OVERRIDES'
    if (.not. read_root()) return
  endif
endif

err_flag = .false.

!-------------------------------------------------------
contains

! read_root reads the top level of root: the lattice, and the tables.

function read_root () result (ok)

logical ok

!

ok = .false.
call root%get_keys (tkeys)
do i = 1, size(tkeys)
  if (any(fel_schema_tables$ == tkeys(i)%key)) then
    call get_value (root, tkeys(i)%key, tab, requested = .false.)
    if (.not. associated(tab)) then
      call out_io (s_error$, r_name, 'IN ' // origin // ', ' // tkeys(i)%key // ' MUST BE A TABLE, [' // &
            tkeys(i)%key // '], AND IS NOT.')
      return
    endif
    if (.not. read_table(tab, tkeys(i)%key)) return
    cycle
  endif
  j = key_index('', tkeys(i)%key)
  if (j == 0) then
    call unknown_key ('', tkeys(i)%key)
    return
  endif
  call fel_schema_set (run, keys(j), root, bad)
  if (bad) return
  if (keys(j)%name == 'lattice') found_lattice = .true.
enddo
ok = .true.

end function read_root

! override_text appends one override to the TOML text and prints it. One that is not
! key=value, or whose key names no table and is not lattice, is refused with the form.
! A bare value for a string or path key, or a string list, is quoted, since the shell
! strips the quotes the user typed.

function override_text (arg) result (ok)

character(*) arg
logical ok
integer ie, id, k
character(:), allocatable :: key, val

!

ok = .false.
ie = index(arg, '=')
if (ie < 2) then
  call out_io (s_error$, r_name, 'AN OVERRIDE ON THE COMMAND LINE IS NOT key=value: ' // trim(arg), &
        'Write table.key=value, as beam.n_particle=4096, or lattice=value.')
  return
endif
key = trim(adjustl(arg(1:ie-1)))
val = trim(adjustl(arg(ie+1:)))
if (val == '') then
  call out_io (s_error$, r_name, 'AN OVERRIDE ON THE COMMAND LINE HAS NO VALUE: ' // trim(arg))
  return
endif
id = index(key, '.', back = .true.)
if (id == 0 .and. key /= 'lattice') then
  call out_io (s_error$, r_name, 'AN OVERRIDE ON THE COMMAND LINE NAMES NO TABLE: ' // trim(arg), &
        'Write table.key=value, as beam.n_particle=4096. Only lattice stands alone.')
  return
endif

if (id == 0) then
  k = key_index('', key)
else
  k = key_index(key(1:id-1), key(id+1:))
endif
if (k > 0) then
  select case (keys(k)%kind)
  case (key_string$, key_path$)
    if (val(1:1) /= '"' .and. val(1:1) /= "'") val = '"' // val // '"'
  case (key_string_list$)
    if (val(1:1) /= '"' .and. val(1:1) /= "'" .and. val(1:1) /= '[') val = '"' // val // '"'
  end select
endif

text = text // key // ' = ' // val // new_line('a')
call out_io (s_info$, r_name, 'Override: ' // key // ' = ' // val)
ok = .true.

end function override_text

! read_table reads every key of one table, and the sub-table [beam.resample].

recursive function read_table (tab, table) result (ok)

type (toml_table), intent(inout) :: tab
character(*) table
logical ok
integer k, jj

!

ok = .false.
call tab%get_keys (ckeys)
do k = 1, size(ckeys)
  if (table == 'beam' .and. ckeys(k)%key == 'resample') then
    call get_value (tab, 'resample', sub, requested = .false.)
    if (.not. associated(sub)) then
      call out_io (s_error$, r_name, 'IN ' // origin // ', [beam] resample MUST BE A TABLE, [beam.resample].')
      return
    endif
    if (.not. read_table(sub, 'beam.resample')) return
    cycle
  endif
  jj = key_index(table, ckeys(k)%key)
  if (jj == 0) then
    call unknown_key (table, ckeys(k)%key)
    return
  endif
  call fel_schema_set (run, keys(jj), tab, bad)
  if (bad) return
enddo
ok = .true.

end function read_table

! key_index finds a key by table and name, 0 when the table has no such key.

function key_index (table, name) result (ix)
character(*) table, name
integer ix, k
ix = 0
do k = 1, size(keys)
  if (keys(k)%table == table .and. keys(k)%name == name) then
    ix = k
    return
  endif
enddo
end function key_index

! unknown_key refuses a key the table does not have, naming the table it belongs in
! when another has it, and else the keys of this table.

subroutine unknown_key (table, name)
character(*) table, name
character(:), allocatable :: list, where
integer k
where = '[' // trim(table) // ']'
if (table == '') where = 'the top level'
do k = 1, size(keys)
  if (keys(k)%name == name) then
    if (keys(k)%table == '') then
      call out_io (s_error$, r_name, 'IN ' // origin // ', ' // name // ' IS NOT A KEY OF ' // where // &
            '. It belongs at the top level.')
    else
      call out_io (s_error$, r_name, 'IN ' // origin // ', ' // name // ' IS NOT A KEY OF ' // where // &
            '. It belongs in [' // trim(keys(k)%table) // '].')
    endif
    return
  endif
enddo
list = ''
do k = 1, size(keys)
  if (keys(k)%table /= table) cycle
  if (list /= '') list = list // ', '
  list = list // trim(keys(k)%name)
enddo
if (table == '') list = list // ', and the tables ' // tables_text()
call out_io (s_error$, r_name, 'IN ' // origin // ', ' // name // ' IS NOT A KEY OF ' // where // &
      '. The keys there are: ' // list // '.')
end subroutine unknown_key

function tables_text () result (txt)
character(:), allocatable :: txt
integer k
txt = ''
do k = 1, size(fel_schema_tables$)
  if (k > 1) txt = txt // ', '
  txt = txt // '[' // trim(fel_schema_tables$(k)) // ']'
enddo
end function tables_text

end subroutine fel_read_input

end module fel_input_mod
