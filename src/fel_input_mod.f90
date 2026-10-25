!+
! Module fel_input_mod
!
! The deck reader of the FEL tracker, kept apart the way Tao keeps its tao_init_*
! files. A deck is one TOML file, read through the schema (fel_schema_mod) into the
! input structures of fel_struct. This module is library code. An embedding program
! may use the reader, or skip it and fill the structures in code.
!
! Every key the deck can set is in the schema's table, and a key the table does not
! name is refused with the keys of its table, so a misspelling cannot pass as a
! default. A key in the wrong table is refused with the table it belongs in. The lattice
! is required, and everything else keeps the declarations' defaults when absent.
!-

module fel_input_mod

use fel_struct
use fel_schema_mod
use tomlf, only: toml_table, toml_key, toml_error, toml_load, get_value

implicit none

contains

!------------------------------------------------------------------------------
!+
! Subroutine fel_read_deck (deck_file, run, err_flag)
!
! Routine to read a TOML deck into the input structures of run. The Bmad switches of the
! [bmad] table set bmad_com and space_charge_com directly, since they are Bmad's own
! globals.
!
! Input:
!   deck_file  -- character(*): The deck.
!
! Output:
!   run        -- fel_run_struct: Input structures filled (global, winit, slicing,
!                   chamber_wake, space_charge, beam_init, resample, bparam, lat_file,
!                   field_file). Nothing else is touched.
!   err_flag   -- logical: Set when the file does not parse, a key is unknown or in the
!                   wrong table, a value has the wrong type or is outside its bounds, or
!                   the lattice is missing.
!-

subroutine fel_read_deck (deck_file, run, err_flag)

type (fel_run_struct), target :: run
type (fel_key_struct), allocatable :: keys(:)
type (toml_table), allocatable :: root
type (toml_table), pointer :: tab, sub
type (toml_key), allocatable :: tkeys(:), ckeys(:)
type (toml_error), allocatable :: terr
character(*) deck_file
logical err_flag, found_lattice, bad
integer i, j, it
character(*), parameter :: r_name = 'fel_read_deck'

!

err_flag = .true.
found_lattice = .false.
call fel_schema_table (keys)

call toml_load (root, trim(deck_file), error = terr)
if (allocated(terr)) then
  call out_io (s_error$, r_name, 'CANNOT READ THE DECK ' // trim(deck_file) // ':', trim(terr%message))
  return
endif

! The top level: the lattice, and the tables.

call root%get_keys (tkeys)
do i = 1, size(tkeys)
  if (any(fel_schema_tables$ == tkeys(i)%key)) then
    call get_value (root, tkeys(i)%key, tab, requested = .false.)
    if (.not. associated(tab)) then
      call out_io (s_error$, r_name, 'IN THE DECK, ' // tkeys(i)%key // ' MUST BE A TABLE, [' // &
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

if (.not. found_lattice) then
  call out_io (s_error$, r_name, 'THE DECK ' // trim(deck_file) // ' NAMES NO LATTICE. ' // &
        'The top-level key lattice is required.')
  return
endif

err_flag = .false.

!-------------------------------------------------------
contains

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
      call out_io (s_error$, r_name, 'IN THE DECK, [beam] resample MUST BE A TABLE, [beam.resample].')
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
      call out_io (s_error$, r_name, 'IN THE DECK, ' // name // ' IS NOT A KEY OF ' // where // &
            '. It belongs at the top level.')
    else
      call out_io (s_error$, r_name, 'IN THE DECK, ' // name // ' IS NOT A KEY OF ' // where // &
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
call out_io (s_error$, r_name, 'IN THE DECK, ' // name // ' IS NOT A KEY OF ' // where // &
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

end subroutine fel_read_deck

end module fel_input_mod
