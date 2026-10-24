!+
! Module fel_h5_mod
!
! Annotated HDF5 writes for the statistics file (fel-physics.md sec-stats).
! No array reaches the file without saying what it is. Every routine here writes a
! dataset and four string attributes with it:
!
!   @unit         The SI or eV unit of the values, or '1' for a dimensionless count.
!   @long_name    Three words, for an axis label or a table heading.
!   @description  One sentence, for a reader who has never seen the file.
!   @axes         The coords/ datasets this dataset's dimensions run over, in h5py
!                 order, comma separated. EG: 'record,slice' means (n_record, n_slice)
!                 with coords/record down and coords/slice across. Every name resolves
!                 to a coords/ dataset, including the trailing label axes of a vector or
!                 a matrix. So a reader never guesses a dimension from its length.
!                 'none' marks a scalar and nothing else. A coordinate names the axis it
!                 defines.
!
! Note: Units are documentation only. The file's numbers are already SI and eV, and a
! reader must not scale by @unit. openPMD's unitSI is the opposite: a factor to apply.
! The two conventions must not be mixed in one file.
!
! The error argument is INTENT(INOUT) and accumulates. A writer calls a run of these
! and tests once at the end. Nothing here stops.
!
! This module also does what HDF5's Fortran layer does not:
!
!   flag datasets    One byte in the file (H5T_STD_I8LE) written from an ordinary
!                    integer buffer. HDF5 has no boolean type, and a flag is never
!                    stored as a float. HDF5 converts on write. Each flag carries
!                    @dtype_hint = 'bool'.
!   string datasets  A fixed-length string array, which h5py reads as an S<n> array.
!                    Bmad's HDF5 layer writes single strings only.
!   text datasets    One string of any length as a scalar dataset (fel_h5_text), for
!                    provenance too big for HDF5's 64 kB attribute cap.
!-

module fel_h5_mod

use hdf5_interface
use h5lt
use hdf5
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
  module procedure fel_h5_flag_rank0
  module procedure fel_h5_flag_rank1
  module procedure fel_h5_flag_rank2
end interface

private annotate

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine annotate (id, name, unit, label, descrip, axes, error)
!
! Routine to attach the four documentation attributes to a dataset just written.
! The dataset is reached through its parent by name.
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

! Never write a zero-length string. HDF5's Fortran wrapper for a string attribute
! writes its terminator into the caller's buffer. A zero-length string walks off the end
! and corrupts the caller's stack. It once overwrote two loop counts in fel_stats_write.
! So a scalar's axes say 'none' rather than nothing.

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
! Routine to write a real dataset with its four documentation attributes.
! Overloaded for ranks 0 through 4 by the interface fel_h5_real.
! Rank 0 writes a true HDF5 scalar, so its shape agrees with @axes = 'none'.
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

! A scalar dataspace and not a one-element array. Bmad's rank-0 writer makes shape (1,).
! That shape contradicts @axes = 'none', and every read would need a [0].

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
integer(hsize_t) dims(4)
integer h5_err
real(rp) val(:,:,:,:)
logical error
character(*) name, unit, label, descrip, axes

! Bmad's dataset writer stops at rank 3, so rank 4 uses the HDF5 Lite API directly.

dims = shape(val)
call H5LTmake_dataset_double_f (id, name, 4, dims, val, h5_err)
error = error .or. (h5_err < 0)
call annotate (id, name, unit, label, descrip, axes, error)

end subroutine fel_h5_real_rank4

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_h5_int_rank<n> (id, name, unit, label, descrip, axes, val, error)
!
! Routine to write an integer dataset with its four documentation attributes.
! Overloaded for ranks 0 through 2 by the interface fel_h5_int.
! Rank 0 is a true HDF5 scalar. Counts are stored as integers and *not* as floats.
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
! and ones. Overloaded for ranks 0 through 2 by the interface fel_h5_flag.
! Rank 0 takes a logical and has no axes argument.
!
! int8 with @unit = '1' is the format's boolean, since HDF5 has no boolean type.
! Each flag also carries @dtype_hint = 'bool'. So the dataset states the convention,
! and a reader does not have to infer it from the root attributes.
!
! The file type is H5T_STD_I8LE and the memory type is the native integer, so HDF5
! converts on write. h5py reads the result as int8.
!
! Input:
!   id       -- integer(hid_t): Group to write into.
!   name     -- character(*): Dataset name.
!   label    -- character(*): Short label, for an axis or a heading.
!   descrip  -- character(*): One-sentence description.
!   axes     -- character(*): Comma separated coords/ names, or ''.
!   val      -- integer: Zeros and ones.
!
! Output:
!   error    -- logical: Accumulates True on any failure.
!-

subroutine fel_h5_flag_rank0 (id, name, label, descrip, val, error)

integer(hid_t) id, s_id, d_id
integer h5_err
integer, target :: buf
logical error, val
character(*) name, label, descrip
type (c_ptr) f_ptr

! A true scalar flag from a logical. The rank-0 case is for input switches, which are
! logicals in the structs. The array cases take 0/1 integer buffers.

call H5Screate_f (H5S_SCALAR_F, s_id, h5_err)
call H5Dcreate_f (id, name, H5T_STD_I8LE, s_id, d_id, h5_err)
if (h5_err < 0) then
  error = .true.
  call H5Sclose_f (s_id, h5_err)
  return
endif
buf = merge(1, 0, val)
f_ptr = c_loc(buf)
call H5Dwrite_f (d_id, H5T_NATIVE_INTEGER, f_ptr, h5_err)
error = error .or. (h5_err < 0)
call H5Dclose_f (d_id, h5_err)
call H5Sclose_f (s_id, h5_err)

call annotate (id, name, '1', label, descrip, '', error)
call hint_bool (id, name, error)

end subroutine fel_h5_flag_rank0

!------------------------------------------------------------------------------

subroutine fel_h5_flag_rank1 (id, name, label, descrip, axes, val, error)

integer(hid_t) id
integer val(:)
logical error
character(*) name, label, descrip, axes

!

call write_flag_dataset (id, name, 1, [size(val)], val, error)
call annotate (id, name, '1', label, descrip, axes, error)
call hint_bool (id, name, error)

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
call hint_bool (id, name, error)

end subroutine fel_h5_flag_rank2

!------------------------------------------------------------------------------
!+
! Subroutine hint_bool (id, name, error)
!
! Routine to mark a flag dataset as this format's boolean. A reader then tests an
! attribute instead of inferring a dtype-and-unit pattern.
!-

subroutine hint_bool (id, name, error)

integer(hid_t) id
integer h5_err
logical error
character(*) name

!

call H5LTset_attribute_string_f (id, name, 'dtype_hint', 'bool', h5_err)
error = error .or. (h5_err < 0)

end subroutine hint_bool

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
! Bmad's HDF5 layer writes single strings only, so the datatype is built here. It is a
! copy of the native character type sized to fel_h5_str_len$. The data go from one
! contiguous c_char buffer, the form h5dwrite_f takes through a c_ptr.
! wavefront_openpmd_mod writes its complex compound the same way.
!
! Input:
!   id       -- integer(hid_t): Group to write into.
!   name     -- character(*): Dataset name.
!   label    -- character(*): Short label, for an axis or a heading.
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
  call H5Screate_f (H5S_SCALAR_F, s_id, h5_err)      ! One string with no axis is a scalar.
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

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_h5_attr_int (obj_id, attrib, val, error)
!
! Routine to attach a true scalar integer attribute to an open object.
!
! Shape expresses arity. One value is a scalar, and a list is an array.
! HDF5's high-level Fortran layer cannot write a scalar attribute.
! H5LTset_attribute_int_f takes a count and an array, so it always makes shape (1,).
! Bmad's hdf5_write_attribute_int_rank0 wraps that call. A reader would then unwrap
! every scalar, and could not tell a one-entry list from a single value.
!
! Input:
!   obj_id   -- integer(hid_t): Open dataset or group to attach to.
!   attrib   -- character(*): Attribute name.
!   val      -- integer: The value.
!
! Output:
!   error    -- logical: Accumulates True on any failure.
!-

subroutine fel_h5_attr_int (obj_id, attrib, val, error)

integer(hid_t) obj_id, s_id, a_id
integer h5_err, val
integer, target :: buf
logical error
character(*) attrib
type (c_ptr) f_ptr

!

call H5Screate_f (H5S_SCALAR_F, s_id, h5_err)
call H5Acreate_f (obj_id, attrib, H5T_NATIVE_INTEGER, s_id, a_id, h5_err)
if (h5_err < 0) then
  error = .true.
  call H5Sclose_f (s_id, h5_err)
  return
endif

buf = val
f_ptr = c_loc(buf)
call H5Awrite_f (a_id, H5T_NATIVE_INTEGER, f_ptr, h5_err)
error = error .or. (h5_err < 0)

call H5Aclose_f (a_id, h5_err)
call H5Sclose_f (s_id, h5_err)

end subroutine fel_h5_attr_int

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_h5_dset_attr_int (id, name, attrib, val, error)
!
! Routine to attach a true scalar integer attribute to a dataset reached by name from
! its parent group. Also see fel_h5_attr_int for why the high-level layer is not used.
!-

subroutine fel_h5_dset_attr_int (id, name, attrib, val, error)

integer(hid_t) id, d_id
integer h5_err, val
logical error
character(*) name, attrib

!

call H5Dopen_f (id, name, d_id, h5_err)
if (h5_err < 0) then
  error = .true.
  return
endif
call fel_h5_attr_int (d_id, attrib, val, error)
call H5Dclose_f (d_id, h5_err)

end subroutine fel_h5_dset_attr_int

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_h5_dset_attr_strs (id, name, attrib, val, error)
!
! Routine to attach a string-array attribute to a dataset reached by name from its
! parent group. Shape expresses arity: a list is an array, length one included.
! Bmad's hdf5_write_attribute_string_rank1 writes the attribute. This routine adds the
! dataset open.
!-

subroutine fel_h5_dset_attr_strs (id, name, attrib, val, error)

integer(hid_t) id, d_id
integer h5_err
logical error, merr
character(*) name, attrib, val(:)

!

call H5Dopen_f (id, name, d_id, h5_err)
if (h5_err < 0) then
  error = .true.
  return
endif
call hdf5_write_attribute_string_rank1 (d_id, attrib, val, merr)
error = error .or. merr
call H5Dclose_f (d_id, h5_err)

end subroutine fel_h5_dset_attr_strs

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_h5_text (id, name, label, descrip, val, error)
!
! Routine to write one string of any length as a scalar string dataset.
!
! The text goes in a dataset and *not* an attribute. HDF5 caps a single attribute at
! 64 kB, and the largest that writes is 65495 bytes (doc/reading-output.md). A scalar
! string dataset takes 3 MB. An echoed namelist or lattice file can pass the cap on an
! ordinary lattice. That failure is silent unless the error is checked.
!
! fel_h5_str writes a fixed-width array for label axes. This routine writes one string
! at its own length. An empty value becomes '(unrecorded)'. A zero-length HDF5 string
! type is invalid, and an empty value would look like a recorded blank.
!
! Input:
!   id       -- integer(hid_t): Group to write into.
!   name     -- character(*): Dataset name.
!   label    -- character(*): Short label, for a heading.
!   descrip  -- character(*): One-sentence description.
!   val      -- character(*): The text. Trailing blanks are dropped.
!
! Output:
!   error    -- logical: Accumulates True on any failure.
!-

subroutine fel_h5_text (id, name, label, descrip, val, error)

integer(hid_t) id, s_id, d_id, t_id
integer h5_err, i, n
logical error
character(*) name, label, descrip, val
character(:), allocatable :: txt
character(kind = c_char), allocatable, target :: buf(:)
type (c_ptr) f_ptr

!

txt = trim(val)
if (txt == '') txt = '(unrecorded)'
n = len(txt)

allocate (buf(n))
do i = 1, n
  buf(i) = txt(i:i)
enddo

call H5Tcopy_f (H5T_NATIVE_CHARACTER, t_id, h5_err)
call H5Tset_size_f (t_id, int(n, size_t), h5_err)
call H5Screate_f (H5S_SCALAR_F, s_id, h5_err)
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

call annotate (id, name, '1', label, descrip, '', error)

end subroutine fel_h5_text

end module fel_h5_mod
