!+
! Module fel_h5_mod
!
! Annotated HDF5 writes for the statistics file (manual sec:stats). ONE RULE: no array
! reaches the file without saying what it is. Every routine here writes a dataset and
! four string attributes with it:
!
!   @unit         The SI or eV unit of the values, or '1' for a dimensionless count.
!   @long_name    Three words, for an axis label or a table heading.
!   @description  One sentence, for a reader who has never seen the file.
!   @axes         The coords/ datasets this dataset's dimensions run over, in h5py
!                 order, comma separated: 'record,s_slice' means (n_record, n_slice)
!                 with coords/record down and coords/s_slice across. EVERY name in it
!                 resolves to a coords/ dataset, including the trailing label axes of a
!                 vector or a matrix, so a reader never guesses a dimension from its
!                 length. 'none' marks a scalar and nothing else: a coordinate names the
!                 axis it defines.
!
! UNITS ARE DOCUMENTATION, NEVER LOAD-BEARING. The file's numbers are already SI and
! eV, so a reader must not scale by @unit. That is the opposite of openPMD's unitSI,
! which is a factor to apply, and the two conventions must not be mixed in one file.
!
! The error argument is INTENT(INOUT) and accumulates: a writer calls a run of these
! and tests once at the end, rather than testing sixty times. Nothing here stops.
!
! What HDF5's Fortran layer does not give us, and this module therefore does:
!
!   flag datasets    One byte in the file (H5T_STD_I8LE) written from an ordinary
!                    integer buffer, since HDF5 has no boolean type and a flag stored
!                    as a float is a lie about what it is. HDF5 converts on write.
!   string datasets  A fixed-length string array, which h5py reads as an S<n> array.
!                    Bmad's HDF5 layer writes single strings only.
!-

module fel_h5_mod

use hdf5_interface
use, intrinsic :: iso_c_binding

implicit none

! Fixed width of a string dataset entry. Bmad element names are 40 characters.

integer, parameter :: fel_h5_str_len$ = 40

interface fel_h5_real
  module procedure fel_h5_real_rank0
  module procedure fel_h5_real_rank1
  module procedure fel_h5_real_rank2
  module procedure fel_h5_real_rank3
  module procedure fel_h5_real_rank4
end interface

interface fel_h5_int
  module procedure fel_h5_int_rank0
  module procedure fel_h5_int_rank1
  module procedure fel_h5_int_rank2
end interface

interface fel_h5_flag
  module procedure fel_h5_flag_rank1
  module procedure fel_h5_flag_rank2
end interface

private annotate

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine annotate (id, name, unit, descrip, axes, error)
!
! Routine to attach the four documentation attributes to a dataset just written.
! The attributes go on the dataset, reached through its parent by name.
!
! Input:
!   id       -- integer(hid_t): Group holding the dataset.
!   name     -- character(*): Dataset name.
!   unit     -- character(*): Unit of the values (documentation only).
!   label    -- character(*): Short label, for an axis or a heading.
!   descrip  -- character(*): One-sentence description.
!   axes     -- character(*): Comma separated coords/ names, or '' for a scalar.
!
! Output:
!   error    -- logical: Accumulates True on any failure.
!-

subroutine annotate (id, name, unit, label, descrip, axes, error)

integer(hid_t) id
integer h5_err
logical error
character(*) name, unit, label, descrip, axes
character(200) ax

!

! NEVER a zero-length string. HDF5's Fortran wrapper for a string attribute writes its
! terminator into the caller's buffer, so a zero-length one walks off the end and
! corrupts the CALLER's stack: it silently overwrote two loop counts in fel_stats_write
! before this was found. An axis or a scalar therefore says 'none' rather than nothing.

ax = axes
if (ax == '') ax = 'none'

call H5LTset_attribute_string_f (id, name, 'unit', unit, h5_err)
error = error .or. (h5_err < 0)
call H5LTset_attribute_string_f (id, name, 'long_name', label, h5_err)
error = error .or. (h5_err < 0)
call H5LTset_attribute_string_f (id, name, 'description', descrip, h5_err)
error = error .or. (h5_err < 0)
call H5LTset_attribute_string_f (id, name, 'axes', trim(ax), h5_err)
error = error .or. (h5_err < 0)

end subroutine annotate

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_h5_real_rank<n> (id, name, unit, label, descrip, axes, val, error)
!
! Routine to write a real dataset with its four documentation attributes. Overloaded for
! ranks 0 through 4 by the interface fel_h5_real. Rank 0 writes a TRUE HDF5 scalar, so a
! scalar's shape and its @axes = 'none' say the same thing.
!
! Input:
!   id       -- integer(hid_t): Group to write into.
!   name     -- character(*): Dataset name.
!   unit     -- character(*): Unit of the values (documentation only).
!   label    -- character(*): Short label, for an axis or a heading.
!   descrip  -- character(*): One-sentence description.
!   axes     -- character(*): Comma separated coords/ names, or ''.
!   val      -- real(rp): The values.
!
! Output:
!   error    -- logical: Accumulates True on any failure.
!-

subroutine fel_h5_real_rank0 (id, name, unit, label, descrip, axes, val, error)

integer(hid_t) id, s_id, d_id
integer h5_err
real(rp) val
real(rp), target :: buf
logical error
character(*) name, unit, label, descrip, axes
type (c_ptr) f_ptr

! A scalar dataspace, not a one-element array: Bmad's rank-0 writer makes shape (1,),
! and then the shape and @axes = 'none' contradict each other and every read needs a [0].

call H5Screate_f (H5S_SCALAR_F, s_id, h5_err)
call H5Dcreate_f (id, name, H5T_NATIVE_DOUBLE, s_id, d_id, h5_err)
if (h5_err < 0) then
  error = .true.
  return
endif
buf = val
f_ptr = c_loc(buf)
call H5Dwrite_f (d_id, H5T_NATIVE_DOUBLE, f_ptr, h5_err)
error = error .or. (h5_err < 0)
call H5Dclose_f (d_id, h5_err)
call H5Sclose_f (s_id, h5_err)

call annotate (id, name, unit, label, descrip, axes, error)

end subroutine fel_h5_real_rank0

!------------------------------------------------------------------------------

subroutine fel_h5_real_rank1 (id, name, unit, label, descrip, axes, val, error)

integer(hid_t) id
real(rp) val(:)
logical error, err
character(*) name, unit, label, descrip, axes

!

call hdf5_write_dataset_real (id, name, val, err)
error = error .or. err
call annotate (id, name, unit, label, descrip, axes, error)

end subroutine fel_h5_real_rank1

!------------------------------------------------------------------------------

subroutine fel_h5_real_rank2 (id, name, unit, label, descrip, axes, val, error)

integer(hid_t) id
real(rp) val(:,:)
logical error, err
character(*) name, unit, label, descrip, axes

!

call hdf5_write_dataset_real (id, name, val, err)
error = error .or. err
call annotate (id, name, unit, label, descrip, axes, error)

end subroutine fel_h5_real_rank2

!------------------------------------------------------------------------------

subroutine fel_h5_real_rank3 (id, name, unit, label, descrip, axes, val, error)

integer(hid_t) id
real(rp) val(:,:,:)
logical error, err
character(*) name, unit, label, descrip, axes

!

call hdf5_write_dataset_real (id, name, val, err)
error = error .or. err
call annotate (id, name, unit, label, descrip, axes, error)

end subroutine fel_h5_real_rank3

!------------------------------------------------------------------------------

subroutine fel_h5_real_rank4 (id, name, unit, label, descrip, axes, val, error)

integer(hid_t) id
real(rp) val(:,:,:,:)
logical error, err
character(*) name, unit, label, descrip, axes

!

call hdf5_write_dataset_real (id, name, val, err)
error = error .or. err
call annotate (id, name, unit, label, descrip, axes, error)

end subroutine fel_h5_real_rank4

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_h5_int_rank<n> (id, name, unit, label, descrip, axes, val, error)
!
! Routine to write an integer dataset with its four documentation attributes. Overloaded
! for ranks 0 through 2 by the interface fel_h5_int, rank 0 being a true HDF5 scalar.
! Counts are integers here, not floats: see the module header.
!
! Input and output as fel_h5_real, with val integer.
!-

subroutine fel_h5_int_rank0 (id, name, unit, label, descrip, axes, val, error)

integer(hid_t) id, s_id, d_id
integer h5_err, val
integer, target :: buf
logical error
character(*) name, unit, label, descrip, axes
type (c_ptr) f_ptr

! A scalar dataspace: see fel_h5_real_rank0.

call H5Screate_f (H5S_SCALAR_F, s_id, h5_err)
call H5Dcreate_f (id, name, H5T_NATIVE_INTEGER, s_id, d_id, h5_err)
if (h5_err < 0) then
  error = .true.
  return
endif
buf = val
f_ptr = c_loc(buf)
call H5Dwrite_f (d_id, H5T_NATIVE_INTEGER, f_ptr, h5_err)
error = error .or. (h5_err < 0)
call H5Dclose_f (d_id, h5_err)
call H5Sclose_f (s_id, h5_err)

call annotate (id, name, unit, label, descrip, axes, error)

end subroutine fel_h5_int_rank0

!------------------------------------------------------------------------------

subroutine fel_h5_int_rank1 (id, name, unit, label, descrip, axes, val, error)

integer(hid_t) id
integer val(:)
logical error, err
character(*) name, unit, label, descrip, axes

!

call hdf5_write_dataset_int (id, name, val, err)
error = error .or. err
call annotate (id, name, unit, label, descrip, axes, error)

end subroutine fel_h5_int_rank1

!------------------------------------------------------------------------------

subroutine fel_h5_int_rank2 (id, name, unit, label, descrip, axes, val, error)

integer(hid_t) id
integer val(:,:)
logical error, err
character(*) name, unit, label, descrip, axes

!

call hdf5_write_dataset_int (id, name, val, err)
error = error .or. err
call annotate (id, name, unit, label, descrip, axes, error)

end subroutine fel_h5_int_rank2

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_h5_flag_rank<n> (id, name, label, descrip, axes, val, error)
!
! Routine to write a one-byte flag dataset from an ordinary integer buffer of zeros
! and ones. Overloaded for ranks 1 and 2 by the interface fel_h5_flag.
!
! int8 with @unit = '1' IS the format's boolean, stated in the file's units_note: HDF5
! has no boolean type, and a flag stored as a float would be a lie about what it is.
!
! The file type is H5T_STD_I8LE and the memory type is the native integer, so HDF5
! converts on write. h5py reads the result as int8, which is the closest thing HDF5
! has to a boolean, and a flag is never stored as a float.
!
! Input:
!   id       -- integer(hid_t): Group to write into.
!   name     -- character(*): Dataset name.
!   descrip  -- character(*): One-sentence description.
!   axes     -- character(*): Comma separated coords/ names, or ''.
!   val      -- integer: Zeros and ones.
!
! Output:
!   error    -- logical: Accumulates True on any failure.
!-

subroutine fel_h5_flag_rank1 (id, name, label, descrip, axes, val, error)

integer(hid_t) id
integer val(:)
logical error
character(*) name, label, descrip, axes

!

call write_flag_dataset (id, name, 1, [size(val)], val, error)
call annotate (id, name, '1', label, descrip, axes, error)

end subroutine fel_h5_flag_rank1

!------------------------------------------------------------------------------

subroutine fel_h5_flag_rank2 (id, name, label, descrip, axes, val, error)

integer(hid_t) id
integer val(:,:)
logical error
character(*) name, label, descrip, axes

!

call write_flag_dataset (id, name, 2, [size(val,1), size(val,2)], val, error)
call annotate (id, name, '1', label, descrip, axes, error)

end subroutine fel_h5_flag_rank2

!------------------------------------------------------------------------------
!+
! Subroutine write_flag_dataset (id, name, rank, dims, val, error)
!
! Routine to create a one-byte integer dataset of the given shape and write val into
! it. val is the caller's array of any rank, taken as a sequence.
!-

subroutine write_flag_dataset (id, name, rank, dims, val, error)

integer(hid_t) id, s_id, d_id
integer rank, dims(:), val(*)
integer(hsize_t) hdims(size(dims))
integer h5_err
logical error
character(*) name

!

hdims = dims

call H5Screate_simple_f (rank, hdims, s_id, h5_err)
if (h5_err < 0) then
  error = .true.
  return
endif

call H5Dcreate_f (id, name, H5T_STD_I8LE, s_id, d_id, h5_err)
if (h5_err < 0) then
  error = .true.
  call H5Sclose_f (s_id, h5_err)
  return
endif

call H5Dwrite_f (d_id, H5T_NATIVE_INTEGER, val(1:product(dims)), hdims, h5_err)
error = error .or. (h5_err < 0)

call H5Dclose_f (d_id, h5_err)
call H5Sclose_f (s_id, h5_err)

end subroutine write_flag_dataset

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_h5_str (id, name, label, descrip, axes, val, error)
!
! Routine to write a fixed-length string dataset, which h5py reads as an S<n> array.
! Bmad's HDF5 layer writes single strings only, so the datatype is built here: a copy
! of the native character type sized to fel_h5_str_len$, written from one contiguous
! c_char buffer (the interoperable form h5dwrite_f takes through a c_ptr, the same
! route wavefront_openpmd_mod uses for its complex compound).
!
! Input:
!   id       -- integer(hid_t): Group to write into.
!   name     -- character(*): Dataset name.
!   descrip  -- character(*): One-sentence description.
!   axes     -- character(*): Comma separated coords/ names, or ''.
!   val(:)   -- character(*): The strings. Longer than fel_h5_str_len$ truncates.
!
! Output:
!   error    -- logical: Accumulates True on any failure.
!-

subroutine fel_h5_str (id, name, label, descrip, axes, val, error)

integer(hid_t) id, s_id, d_id, t_id
integer h5_err, i, j, n
integer(hsize_t) hdims(1)
logical error
character(*) name, label, descrip, axes, val(:)
character(kind = c_char), allocatable, target :: buf(:)
type (c_ptr) f_ptr

!

n = size(val)
hdims = n
allocate (buf(fel_h5_str_len$ * max(n, 1)))
buf = c_char_' '

do i = 1, n
  do j = 1, min(len_trim(val(i)), fel_h5_str_len$)
    buf((i-1) * fel_h5_str_len$ + j) = val(i)(j:j)
  enddo
enddo

call H5Tcopy_f (H5T_NATIVE_CHARACTER, t_id, h5_err)
call H5Tset_size_f (t_id, int(fel_h5_str_len$, size_t), h5_err)
if (n == 1 .and. axes == '') then
  call H5Screate_f (H5S_SCALAR_F, s_id, h5_err)      ! One string with no axis IS a scalar.
else
  call H5Screate_simple_f (1, hdims, s_id, h5_err)
endif
call H5Dcreate_f (id, name, t_id, s_id, d_id, h5_err)
if (h5_err < 0) then
  error = .true.
  return
endif

f_ptr = c_loc(buf(1))
call H5Dwrite_f (d_id, t_id, f_ptr, h5_err)
error = error .or. (h5_err < 0)

call H5Dclose_f (d_id, h5_err)
call H5Sclose_f (s_id, h5_err)
call H5Tclose_f (t_id, h5_err)

call annotate (id, name, '1', label, descrip, axes, error)

end subroutine fel_h5_str

end module fel_h5_mod
