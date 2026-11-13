!+
! Module fel_print_mod
!
! The layout of what a run prints: the logo, and boxes of rows in one or two columns
! with titled rules between them. A row is composed in single-byte characters at fixed
! widths, with stand-ins for the borders and for the two characters a unit or the logo
! may need, the micro sign and the middle dot. fel_print_line then writes the row with
! the box-drawing characters of Unicode or with ASCII. A box-drawing character is three
! bytes and one column, so no width is ever counted on the text that is written.
!
! Unicode is used where the locale names UTF-8 (LC_ALL, LC_CTYPE, LANG, in that order),
! and the program's --ascii and --unicode set it outright. The two forms have the same
! columns.
!-

module fel_print_mod

use sim_utils

implicit none

! The box: fel_box_width$ columns from the left border to the right one, with the middle
! border of a two-column box at fel_box_half$. Every line starts with one blank, so a box
! ends in column 80.

integer, parameter :: fel_box_width$ = 78
integer, parameter :: fel_box_half$ = 39

! Stand-ins, single bytes that no key, value or file name holds.

character(1), parameter :: fel_bar$ = achar(1)        ! A vertical border.
character(1), parameter :: fel_micro$ = achar(2)      ! The micro sign of a unit.
character(1), parameter :: fel_dot$ = achar(3)        ! The middle dot of the logo.
character(1), parameter :: fel_subrule$ = achar(4)    ! First character of a column entry that is a rule.

! The kinds of rule.

integer, parameter :: fel_rule_top$ = 1, fel_rule_mid$ = 2, fel_rule_bottom$ = 3

!+
! Struct fel_print_struct
!
! How the run prints. unicode selects the box-drawing characters, verbose prints the
! setup's information lines on the screen as well as in the log, tty is true where
! standard output is a terminal, and progress is how the tracking reports: bar, one line
! rewritten in place, which a terminal takes unless told otherwise, table, a row at each
! element end and every two seconds inside a long element, or none.
!-

type fel_print_struct
  logical :: unicode = .false.
  logical :: verbose = .false.
  logical :: tty = .false.
  character(8) :: progress = 'table'
end type

type (fel_print_struct), save :: fel_print_com

!+
! Struct fel_column_struct
!
! One column of a two-column box, or the rows of a one-column box, gathered before any
! is written so that the keys can be padded to the longest. An entry is a key, a value
! and a provenance mark, or a rule across the column when its key starts with
! fel_subrule$ and its value is the rule's title.
!-

type fel_column_struct
  character(40), allocatable :: key(:)
  character(120), allocatable :: val(:)
  character(1), allocatable :: mark(:)
  integer :: n = 0
end type

private :: utf8, rule_cells, locale_is_utf8

interface
  function c_isatty (fd) bind(c, name = 'isatty') result (is)
    use, intrinsic :: iso_c_binding, only: c_int
    integer(c_int), value :: fd
    integer(c_int) is
  end function
  function c_gethostname (name, n) bind(c, name = 'gethostname') result (ierr)
    use, intrinsic :: iso_c_binding, only: c_int, c_char, c_size_t
    character(kind = c_char) name(*)
    integer(c_size_t), value :: n
    integer(c_int) ierr
  end function
end interface

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_print_setup (force_ascii, force_unicode, verbose, progress)
!
! Routine to choose how the run prints. Unicode is taken where the locale names UTF-8,
! and either argument sets it outright, ASCII winning where both are given.
!
! Input:
!   force_ascii   -- logical: The command line asked for ASCII.
!   force_unicode -- logical: The command line asked for Unicode.
!   verbose       -- logical: Print the setup's information lines on the screen.
!   progress      -- character(*): bar, table or none, or blank for the bar on a terminal
!                      and the table elsewhere.
!-

subroutine fel_print_setup (force_ascii, force_unicode, verbose, progress)

use, intrinsic :: iso_c_binding, only: c_int

logical force_ascii, force_unicode, verbose
character(*) progress

!

fel_print_com%unicode = locale_is_utf8()
if (force_unicode) fel_print_com%unicode = .true.
if (force_ascii) fel_print_com%unicode = .false.
fel_print_com%verbose = verbose
fel_print_com%tty = (c_isatty(1_c_int) /= 0)
fel_print_com%progress = progress
if (progress == '') fel_print_com%progress = merge('bar  ', 'table', fel_print_com%tty)

end subroutine fel_print_setup

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function locale_is_utf8 () result (is_utf8)
!
! Routine to say whether the locale names UTF-8: the first of LC_ALL, LC_CTYPE and LANG
! that is set decides, as the C library reads them.
!
! Output:
!   is_utf8 -- logical: True where that variable holds UTF-8 or utf8 in any case.
!-

function locale_is_utf8 () result (is_utf8)

logical is_utf8
character(64) val
character(8), parameter :: names(3) = [character(8) :: 'LC_ALL', 'LC_CTYPE', 'LANG']
integer i, n, stat

!

is_utf8 = .false.
do i = 1, size(names)
  call get_environment_variable (trim(names(i)), val, n, stat)
  if (stat /= 0 .or. n == 0) cycle
  call str_upcase (val, val)
  is_utf8 = (index(val, 'UTF-8') > 0 .or. index(val, 'UTF8') > 0)
  return
enddo

end function locale_is_utf8

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_print_host () result (host)
!
! Routine to return the name of the machine the run is on, which the logo states. A batch
! system's logs from many nodes are told apart by it.
!
! Output:
!   host -- character(64): The host name, blank where the system gives none.
!-

function fel_print_host () result (host)

use, intrinsic :: iso_c_binding, only: c_char, c_size_t, c_null_char

character(64) host
character(kind = c_char) buf(65)
integer i

!

host = ''
buf = c_null_char
if (c_gethostname(buf, int(size(buf) - 1, c_size_t)) /= 0) return
do i = 1, len(host)
  if (buf(i) == c_null_char) exit
  host(i:i) = buf(i)
enddo

end function fel_print_host

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function utf8 (name) result (c)
!
! Routine to return a box-drawing character by name as its UTF-8 bytes, or its ASCII
! stand-in where the run prints ASCII.
!
! Input:
!   name -- character(*): h and v the lines, tl tr bl br the rounded corners, lt and rt
!             the left and right joints, tt and bt the joints that open down and up, x
!             the cross, tick_down and tick_up the half lines of the logo.
!
! Output:
!   c    -- character(:), allocatable: The character's bytes.
!-

function utf8 (name) result (c)

character(*) name
character(:), allocatable :: c

!

if (.not. fel_print_com%unicode) then
  select case (name)
  case ('h');                  c = '-'
  case ('v', 'tick_down', 'tick_up');  c = '|'
  case default;                c = '+'
  end select
  return
endif

select case (name)
case ('h');          c = achar(226) // achar(148) // achar(128)
case ('v');          c = achar(226) // achar(148) // achar(130)
case ('tl');         c = achar(226) // achar(149) // achar(173)
case ('tr');         c = achar(226) // achar(149) // achar(174)
case ('bl');         c = achar(226) // achar(149) // achar(176)
case ('br');         c = achar(226) // achar(149) // achar(175)
case ('lt');         c = achar(226) // achar(148) // achar(156)
case ('rt');         c = achar(226) // achar(148) // achar(164)
case ('tt');         c = achar(226) // achar(148) // achar(172)
case ('bt');         c = achar(226) // achar(148) // achar(180)
case ('x');          c = achar(226) // achar(148) // achar(188)
case ('tick_down');  c = achar(226) // achar(149) // achar(183)
case ('tick_up');    c = achar(226) // achar(149) // achar(181)
case default;        c = '?'
end select

end function utf8

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_print_line (line)
!
! Routine to write one composed line: the stand-ins become the characters of the form
! the run prints in, and the line goes out through out_io, to the screen and the log.
!
! Input:
!   line -- character(*): The line in single-byte characters, with fel_bar$ for a
!             border, fel_micro$ for the micro sign and fel_dot$ for the middle dot.
!-

subroutine fel_print_line (line)

character(*) line
character(:), allocatable :: out
integer i
character(*), parameter :: r_name = 'lucifer'

!

out = ''
do i = 1, len_trim(line)
  select case (line(i:i))
  case (fel_bar$)
    out = out // utf8('v')
  case (fel_micro$)
    if (fel_print_com%unicode) then
      out = out // achar(194) // achar(181)
    else
      out = out // 'u'
    endif
  case (fel_dot$)
    if (fel_print_com%unicode) then
      out = out // achar(194) // achar(183)
    else
      out = out // '|'
    endif
  case default
    out = out // line(i:i)
  end select
enddo

call out_io (s_blank$, r_name, out)

end subroutine fel_print_line

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function rule_cells (title, i_first, i_last, left) result (out)
!
! Routine to compose the part of a rule between two columns of the box: the character at
! its left end, the title after two lines and a blank, and lines to its right end.
!
! Input:
!   title   -- character(*): The title, blank for none.
!   i_first -- integer: Column of the left end, counted from the box's left border as 0.
!   i_last  -- integer: Column of the last line of this part.
!   left    -- character(*): The name of the character at the left end.
!
! Output:
!   out     -- character(:), allocatable: The bytes of the part.
!-

function rule_cells (title, i_first, i_last, left) result (out)

character(*) title, left
integer i_first, i_last, n_used, i
character(:), allocatable :: out

!

out = utf8(left)
n_used = 1
if (title /= '') then
  out = out // utf8('h') // utf8('h') // ' ' // trim(title) // ' '
  n_used = n_used + 4 + len_trim(title)
endif
do i = i_first + n_used, i_last
  out = out // utf8('h')
enddo

end function rule_cells

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_box_rule (kind, title, title2, above, below)
!
! Routine to write a rule of a box: its top, a rule between two parts of it, or its
! bottom. The middle border joins the rule where a two-column part is above it, below it
! or both.
!
! Input:
!   kind   -- integer: fel_rule_top$, fel_rule_mid$ or fel_rule_bottom$.
!   title  -- character(*): Title of the part below, at the left. Blank for none.
!   title2 -- character(*): Title at the middle border, of the right column below. Blank
!               for none.
!   above  -- logical: The part above has two columns.
!   below  -- logical: The part below has two columns.
!-

subroutine fel_box_rule (kind, title, title2, above, below)

integer kind
character(*) title, title2
logical above, below
character(:), allocatable :: out, mid
character(2) left, right
character(*), parameter :: r_name = 'lucifer'

!

select case (kind)
case (fel_rule_top$);     left = 'tl';  right = 'tr'
case (fel_rule_bottom$);  left = 'bl';  right = 'br'
case default;             left = 'lt';  right = 'rt'
end select

if (above .and. below) then
  mid = 'x'
elseif (below) then
  mid = 'tt'
elseif (above) then
  mid = 'bt'
else
  mid = 'h'
endif

if (above .or. below) then
  out = ' ' // rule_cells(title, 0, fel_box_half$ - 1, trim(left)) // &
        rule_cells(title2, fel_box_half$, fel_box_width$ - 1, mid) // utf8(trim(right))
else
  out = ' ' // rule_cells(title, 0, fel_box_width$ - 1, trim(left)) // utf8(trim(right))
endif

call out_io (s_blank$, r_name, out)

end subroutine fel_box_rule

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_box_full (text)
!
! Routine to write one row across the whole box. A text longer than the box is written
! whole and passes the right border, which a long file name may do.
!
! Input:
!   text -- character(*): The row's text.
!-

subroutine fel_box_full (text)

character(*) text
character(fel_box_width$ - 2) body

!

if (len_trim(text) > len(body)) then
  call fel_print_line (' ' // fel_bar$ // ' ' // trim(text))
else
  body = text
  call fel_print_line (' ' // fel_bar$ // ' ' // body // fel_bar$)
endif

end subroutine fel_box_full

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_column_add (col, key, val, mark)
!
! Routine to add a row to a column: a key, its value and its provenance mark.
!
! Input:
!   col  -- fel_column_struct: The column.
!   key  -- character(*): The key, or the name of a quantity.
!   val  -- character(*): The value with its unit.
!   mark -- character(*), optional: The provenance mark, one character. Blank by default.
!
! Output:
!   col  -- fel_column_struct: With the row added.
!-

subroutine fel_column_add (col, key, val, mark)

type (fel_column_struct) col
character(*) key, val
character(*), optional :: mark
character(40), allocatable :: k(:)
character(120), allocatable :: v(:)
character(1), allocatable :: m(:)

!

if (.not. allocated(col%key)) allocate (col%key(16), col%val(16), col%mark(16))
if (col%n == size(col%key)) then
  allocate (k(2 * col%n), v(2 * col%n), m(2 * col%n))
  k(1:col%n) = col%key;  v(1:col%n) = col%val;  m(1:col%n) = col%mark
  call move_alloc (k, col%key);  call move_alloc (v, col%val);  call move_alloc (m, col%mark)
endif

col%n = col%n + 1
col%key(col%n) = key
col%val(col%n) = val
col%mark(col%n) = ' '
if (present(mark)) col%mark(col%n) = mark

end subroutine fel_column_add

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_column_rule (col, title)
!
! Routine to add a titled rule across a column, which starts a box of its own within the
! right column of a two-column part.
!
! Input:
!   col   -- fel_column_struct: The column.
!   title -- character(*): The title.
!
! Output:
!   col   -- fel_column_struct: With the rule added.
!-

subroutine fel_column_rule (col, title)

type (fel_column_struct) col
character(*) title

!

call fel_column_add (col, fel_subrule$, title)

end subroutine fel_column_rule

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_column_text (col, i, width, key_width) result (text)
!
! Routine to compose row i of a column at a width: the key padded to key_width, the
! value, and the mark in the last column but one. A value too long for its place is
! written whole and takes the mark's place, and the text is cut at the width.
!
! Input:
!   col       -- fel_column_struct: The column.
!   i         -- integer: The row. Past the column's last, the text is blank.
!   width     -- integer: Columns of the text.
!   key_width -- integer: Columns the key is padded to.
!
! Output:
!   text      -- character(width): The row.
!-

function fel_column_text (col, i, width, key_width) result (text)

type (fel_column_struct) col
integer i, width, key_width, n_val
character(width) text
character(:), allocatable :: row

!

text = ''
if (i > col%n) return

row = col%key(i)(1:max(key_width, len_trim(col%key(i)))) // ' ' // trim(col%val(i))
text = row
if (col%mark(i) /= ' ') then
  n_val = len(row)
  if (n_val <= width - 3) then
    text(width-1:width-1) = col%mark(i)
  else
    text = row // ' ' // col%mark(i)
  endif
endif

end function fel_column_text

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_column_key_width (col) result (width)
!
! Routine to return the columns the keys of a column are padded to: the longest key of
! the column, the rules aside.
!
! Input:
!   col   -- fel_column_struct: The column.
!
! Output:
!   width -- integer: The longest key's length, at least 1.
!-

function fel_column_key_width (col) result (width)

type (fel_column_struct) col
integer width, i

!

width = 1
do i = 1, col%n
  if (col%key(i)(1:1) == fel_subrule$) cycle
  width = max(width, len_trim(col%key(i)))
enddo

end function fel_column_key_width

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_box_columns (left, right, key_width)
!
! Routine to write a two-column part of a box from its two columns, row beside row, the
! shorter column ending in blank rows. A rule in the right column is written from the
! middle border to the right one.
!
! Input:
!   left      -- fel_column_struct: The left column.
!   right     -- fel_column_struct: The right column.
!   key_width -- integer, optional: Columns the keys of both are padded to. By default
!                  each column pads to its own longest key.
!-

subroutine fel_box_columns (left, right, key_width)

type (fel_column_struct) left, right
integer, optional :: key_width
integer i, kw_l, kw_r, w_l, w_r
character(:), allocatable :: out
character(*), parameter :: r_name = 'lucifer'

!

kw_l = fel_column_key_width(left)
kw_r = fel_column_key_width(right)
if (present(key_width)) then
  kw_l = key_width;  kw_r = key_width
endif
w_l = fel_box_half$ - 2
w_r = fel_box_width$ - fel_box_half$ - 2

do i = 1, max(left%n, right%n)
  if (i <= right%n) then
    if (right%key(i)(1:1) == fel_subrule$) then

      ! The left row with its border, then the rule in the box-drawing characters. This
      ! line is written here and not by fel_print_line, which would count the rule's
      ! bytes as columns nowhere, but the left text may hold a stand-in.

      out = ' ' // utf8('v') // ' ' // stand_ins(fel_column_text(left, i, w_l, kw_l)) // &
            rule_cells(trim(right%val(i)), fel_box_half$, fel_box_width$ - 1, 'lt') // utf8('rt')
      call out_io (s_blank$, r_name, out)
      cycle
    endif
  endif
  call fel_print_line (' ' // fel_bar$ // ' ' // fel_column_text(left, i, w_l, kw_l) // fel_bar$ // ' ' // &
                       fel_column_text(right, i, w_r, kw_r) // fel_bar$)
enddo

end subroutine fel_box_columns

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_si_micro (str) result (out)
!
! Routine to put the stand-in fel_micro$ where a number with its unit, as fel_si_str
! writes it, carries the prefix u, so that fel_print_line writes the micro sign in
! Unicode and u in ASCII.
!
! Input:
!   str -- character(*): The number and its unit, such as "48.824 uJ".
!
! Output:
!   out -- character(len(str)): The same text with the stand-in for the prefix.
!-

function fel_si_micro (str) result (out)

character(*) str
character(len(str)) out
integer i

!

out = str
i = len_trim(out)
if (i < 3) return
if (out(i-2:i-1) == ' u') out(i-1:i-1) = fel_micro$

end function fel_si_micro

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function stand_ins (text) result (out)
!
! Routine to replace the micro sign's stand-in in a text by the character of the form
! the run prints in, keeping the text's trailing blanks.
!
! Input:
!   text -- character(*): The text, at its full width.
!
! Output:
!   out  -- character(:), allocatable: The text as written.
!-

function stand_ins (text) result (out)

character(*) text
character(:), allocatable :: out
integer i

!

out = ''
do i = 1, len(text)
  if (text(i:i) == fel_micro$) then
    if (fel_print_com%unicode) then
      out = out // achar(194) // achar(181)
    else
      out = out // 'u'
    endif
  else
    out = out // text(i:i)
  endif
enddo

end function stand_ins

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_print_logo (line1, line2, line3)
!
! Routine to write the logo with three lines of text to its right: the name in
! line-drawn letters with rounded corners in Unicode, and in capitals in ASCII.
!
! Input:
!   line1, line2, line3 -- character(*): The text beside the three rows. fel_dot$
!                            stands for the middle dot.
!-

subroutine fel_print_logo (line1, line2, line3)

character(*) line1, line2, line3
character(*), parameter :: r_name = 'lucifer'

! The letters, a character a column: r 7 L j the corners, - and | the lines, comma and
! apostrophe the half lines that open down and up, + the left joint, T the joint that
! opens down.

character(39), parameter :: rows(3) = [character(39) :: &
      "  ,     ,  ,  r---  ,  r---  r---  r--7", &
      "  |     |  |  |     |  +--   +--   +-Tj", &
      "  L---  L--j  L---  '  '     L---  ' L-"]
character(120) text(3)
character(:), allocatable :: out
integer ir, i

!

text(1) = line1;  text(2) = line2;  text(3) = line3

if (.not. fel_print_com%unicode) then
  call fel_print_line ('  L U C I F E R   ' // trim(text(1)))
  call fel_print_line ('                  ' // trim(text(2)))
  call fel_print_line ('                  ' // trim(text(3)))
  return
endif

do ir = 1, 3
  out = ''
  do i = 1, len(rows(ir))
    select case (rows(ir)(i:i))
    case ('r');  out = out // utf8('tl')
    case ('7');  out = out // utf8('tr')
    case ('L');  out = out // utf8('bl')
    case ('j');  out = out // utf8('br')
    case ('-');  out = out // utf8('h')
    case ('|');  out = out // utf8('v')
    case (',');  out = out // utf8('tick_down')
    case ("'");  out = out // utf8('tick_up')
    case ('+');  out = out // utf8('lt')
    case ('T');  out = out // utf8('tt')
    case default;  out = out // rows(ir)(i:i)
    end select
  enddo
  out = out // '   '
  do i = 1, len_trim(text(ir))
    if (text(ir)(i:i) == fel_dot$) then
      out = out // achar(194) // achar(183)
    else
      out = out // text(ir)(i:i)
    endif
  enddo
  call out_io (s_blank$, r_name, out)
enddo

end subroutine fel_print_logo

end module
