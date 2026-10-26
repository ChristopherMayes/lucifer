!+
! Module wavefront_openpmd_mod
!
! openPMD EXT_Wavefront I/O for wavefront_struct: the openPMD-standard base plus the
! Wavefront extension (openPMD-standard, branch upcoming-2.0.0, EXT_Wavefront.md).
! The standard document is authoritative for this format. This module and the
! benchmark's h5py reader are each checked against its text.
!
! Layout (also see fel-physics.md sec-field-set):
!
!   /                        openPMD "2.0.0", openPMDextension "Wavefront",
!                            basePath "/data/%T/", meshesPath "meshes/",
!                            iterationEncoding "groupBased", iterationFormat "/data/%T/"
!   /data/1/                 One iteration per file. time = 0, dt = 0, timeUnitSI = 1.
!   /data/1/meshes/electricField
!                            The one mesh record. The extension's heading puts its
!                            attributes on the mesh record, and its body says "series".
!                            The record is used, since photonEnergy is a property of one
!                            field. The attributes:
!                            geometry 'cartesian', axisLabels ['z','y','x'],
!                            gridSpacing [dz,dy,dx], gridGlobalOffset [z0,-gmax,-gmax]
!                            (z0 is 0 for a whole window and the first slice's z for a
!                            frame cut to a slice range),
!                            gridUnitSI [1,1,1], gridUnitDimension (a length per axis),
!                            unitDimension (1,1,-3,-1,0,0,0) (V/m),
!                            timeOffset 0, photonEnergy [J], temporalDomain 'time',
!                            spatialDomain 'r', zCoordinate [m].
!   .../electricField/x      Complex compound {r,i} dataset, which h5py reads natively.
!                            Component attributes: unitSI 1 and position [0,0,0], the
!                            sample's offset within its cell in units of gridSpacing.
!   .../electricField/y      Present only when the wavefront carries Ey. Both transverse
!                            polarizations are components in one file. The Genesis
!                            format writes one per file. The z component is never
!                            written: a paraxial code has none, and an absent component
!                            is ordinary openPMD.
!
! The dataset is stored exactly as the Fortran (nx, ny, nslice) array. The HDF5 Fortran
! API records that as a C-order (nslice, ny, nx) dataspace. This needs no copy and gives
! numpy per-slice access directly. axisLabels declare that stored order. The slice axis
! is the one labeled z (third in the logical x,y,z reading). Slices are simultaneous, so
! they are a mesh axis and never the openPMD iteration.
!
! Reading implements exactly what writing produces and refuses the rest:
! temporalDomain 'frequency', spatialDomain 'k', an axisLabels order other than
! (z,y,x), and any missing required attribute.
!
! wavefront_write_openpmd writes one wavefront whole. A file too large to hold in memory,
! or one whose slice count is not known until the run ends, is written a few planes at a
! time: wavefront_openpmd_create makes the file with its attributes and empty component
! datasets, wavefront_openpmd_write_planes fills a run of planes along z, and
! wavefront_openpmd_read_planes reads one back. The escaped-field bank and the pulse
! (fel_io_mod) are written this way.
!-

module wavefront_openpmd_mod

use wavefront_mod
use hdf5_openpmd_mod
use sim_utils

implicit none

character(*), parameter, private :: iter_path = 'data/1'
character(*), parameter, private :: mesh_path = 'data/1/meshes/electricField'

private write_mesh_attributes, write_component_attributes

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function file_is_openpmd (file_name) result (is_pmd)
!
! Routine to test whether a file is an openPMD file. The import path refuses a file on
! this test. It holds for any openPMD file, particles as well as fields. An openPMD file
! carries the required root attribute "openPMD", and a Genesis dump of either kind has
! none.
!
! Input:
!   file_name   -- character(*): File to probe.
!
! Output:
!   is_pmd      -- logical: True if the file carries the openPMD root attribute.
!-

function file_is_openpmd (file_name) result (is_pmd)

character(*) file_name
logical is_pmd, err
integer(hid_t) f_id
integer h5_err
character(40) version

!

is_pmd = .false.
call hdf5_open_file (file_name, 'READ', f_id, err, .false.)
if (err) return
call hdf5_read_attribute_string (f_id, 'openPMD', version, err, .false.)
is_pmd = .not. err
call h5fclose_f (f_id, h5_err)

end function file_is_openpmd

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_write_openpmd (wf, file_name, err_flag, z_offset)
!
! Routine to write one wavefront as an openPMD EXT_Wavefront file, in the layout of the
! module header. The field is already in V/m (temporalDomain 'time'), so unitSI = 1
! throughout. zCoordinate is wf%s_position, the plane's beamline position. It is taken
! from the wavefront so the file cannot disagree with the structure it came from.
!
! Input:
!   wf         -- wavefront_struct: The field. Ey written when allocated.
!   file_name  -- character(*): Output file.
!   z_offset   -- real(rp), optional: Where the first record sits along the window [m],
!                   the z of gridGlobalOffset. A field cut to a slice range passes the
!                   first slice's z. The mesh then lands where the whole window would
!                   have put those slices. Default 0, the whole window.
!
! Output:
!   err_flag   -- logical: Set True on error, False otherwise.
!-

subroutine wavefront_write_openpmd (wf, file_name, err_flag, z_offset)

type (wavefront_struct), target :: wf
character(*) file_name
real(rp), optional :: z_offset
integer(hid_t) f_id, it_id, m_id, complex_t
integer h5_err
logical err_flag, err
character(*), parameter :: r_name = 'wavefront_write_openpmd'

!

err_flag = .true.

call hdf5_open_file (file_name, 'WRITE', f_id, err);  if (err) return

! Root: the openPMD base standard's required series attributes.

call hdf5_write_attribute_string (f_id, 'openPMD', '2.0.0', err)
call hdf5_write_attribute_string (f_id, 'openPMDextension', 'Wavefront', err)
call hdf5_write_attribute_string (f_id, 'basePath', '/data/%T/', err)
call hdf5_write_attribute_string (f_id, 'meshesPath', 'meshes/', err)
call hdf5_write_attribute_string (f_id, 'iterationEncoding', 'groupBased', err)
call hdf5_write_attribute_string (f_id, 'iterationFormat', '/data/%T/', err)

! One iteration. Slices are a mesh axis, so the iteration carries no time structure.

call h5gcreate_f (f_id, 'data', it_id, h5_err);  call h5gclose_f (it_id, h5_err)
call h5gcreate_f (f_id, trim(iter_path), it_id, h5_err)
call hdf5_write_attribute_real (it_id, 'time', 0.0_rp, err)
call hdf5_write_attribute_real (it_id, 'dt', 0.0_rp, err)
call hdf5_write_attribute_real (it_id, 'timeUnitSI', 1.0_rp, err)
call h5gclose_f (it_id, h5_err)
call h5gcreate_f (f_id, 'data/1/meshes', it_id, h5_err);  call h5gclose_f (it_id, h5_err)

! The mesh record, with the base standard's mesh attributes and the Wavefront
! extension's required set on the record (module header).

call h5gcreate_f (f_id, trim(mesh_path), m_id, h5_err)
if (h5_err < 0) then
  call out_io (s_error$, r_name, 'CANNOT CREATE MESH GROUP IN ' // trim(file_name))
  return
endif

call write_mesh_attributes (m_id, wf, real_option(0.0_rp, z_offset), wf%s_position, err)
if (err) return

call pmd_init_compound_complex (complex_t)

call pmd_write_complex_to_dataset (m_id, 'x', complex_t, 'x', unit_V_per_m, wf%Ex, err)
if (err) then
  call pmd_kill_compound_complex (complex_t)
  return
endif
call write_component_attributes (m_id, 'x', err)

if (allocated(wf%Ey)) then
  call pmd_write_complex_to_dataset (m_id, 'y', complex_t, 'y', unit_V_per_m, wf%Ey, err)
  if (err) then
    call pmd_kill_compound_complex (complex_t)
    return
  endif
  call write_component_attributes (m_id, 'y', err)
endif

call pmd_kill_compound_complex (complex_t)
call h5gclose_f (m_id, h5_err)
call h5fclose_f (f_id, h5_err)

err_flag = .false.

end subroutine wavefront_write_openpmd

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_openpmd_create (file_name, wf, nz, f_id, m_id, err_flag, z_offset, extendable, z_coordinate)
!
! Routine to create an openPMD EXT_Wavefront file in the layout of the module header,
! with its attributes taken from wf and its component datasets created empty, for a
! caller that fills the planes with wavefront_openpmd_write_planes. The x component is
! always created, and y when wf carries Ey. The caller closes m_id and f_id.
!
! The grid and the wavelength come from wf, and nz is the slice count of the file, which
! need not be wf's. With extendable True the datasets are chunked one plane at a time
! and grow without bound along z, so nz is their first extent and a plane written past
! it extends them. Otherwise the datasets are contiguous at nz planes, as
! wavefront_write_openpmd writes them.
!
! Input:
!   file_name     -- character(*): Output file.
!   wf            -- wavefront_struct: The grid, the wavelength and the components.
!   nz            -- integer: Slice count of the datasets.
!   z_offset      -- real(rp), optional: Where the first plane sits along the pulse [m],
!                      the z of gridGlobalOffset. Default 0.
!   extendable    -- logical, optional: Chunked datasets that grow along z? Default False.
!   z_coordinate  -- real(rp), optional: The beamline position written as zCoordinate [m].
!                      Default is wf%s_position.
!
! Output:
!   f_id          -- integer(hid_t): The open file.
!   m_id          -- integer(hid_t): The open mesh record group, the parent of x and y.
!   err_flag      -- logical: Set True on error, False otherwise.
!-

subroutine wavefront_openpmd_create (file_name, wf, nz, f_id, m_id, err_flag, z_offset, &
                                     extendable, z_coordinate)

type (wavefront_struct), target :: wf
character(*) file_name
real(rp), optional :: z_offset, z_coordinate
logical, optional :: extendable
integer nz
integer(hid_t) f_id, m_id, it_id, complex_t, s_id, d_id, p_id
integer(hsize_t) dims(3), maxdims(3)
integer h5_err
logical err_flag, err, grow
character(*), parameter :: r_name = 'wavefront_openpmd_create'

!

err_flag = .true.
grow = logic_option(.false., extendable)

call hdf5_open_file (file_name, 'WRITE', f_id, err);  if (err) return

call hdf5_write_attribute_string (f_id, 'openPMD', '2.0.0', err)
call hdf5_write_attribute_string (f_id, 'openPMDextension', 'Wavefront', err)
call hdf5_write_attribute_string (f_id, 'basePath', '/data/%T/', err)
call hdf5_write_attribute_string (f_id, 'meshesPath', 'meshes/', err)
call hdf5_write_attribute_string (f_id, 'iterationEncoding', 'groupBased', err)
call hdf5_write_attribute_string (f_id, 'iterationFormat', '/data/%T/', err)

call h5gcreate_f (f_id, 'data', it_id, h5_err);  call h5gclose_f (it_id, h5_err)
call h5gcreate_f (f_id, trim(iter_path), it_id, h5_err)
call hdf5_write_attribute_real (it_id, 'time', 0.0_rp, err)
call hdf5_write_attribute_real (it_id, 'dt', 0.0_rp, err)
call hdf5_write_attribute_real (it_id, 'timeUnitSI', 1.0_rp, err)
call h5gclose_f (it_id, h5_err)
call h5gcreate_f (f_id, 'data/1/meshes', it_id, h5_err);  call h5gclose_f (it_id, h5_err)

call h5gcreate_f (f_id, trim(mesh_path), m_id, h5_err)
if (h5_err < 0) then
  call out_io (s_error$, r_name, 'CANNOT CREATE MESH GROUP IN ' // trim(file_name))
  return
endif

call write_mesh_attributes (m_id, wf, real_option(0.0_rp, z_offset), &
                            real_option(wf%s_position, z_coordinate), err)
if (err) return

! The component datasets, empty. The Fortran (nx, ny, nz) shape is what the HDF5 Fortran
! API records as a C-order (nz, ny, nx) dataspace, as the whole-file writer's is. An
! extendable dataset must be chunked, and a chunk of one plane is the unit the planes
! are written in.

dims = [size(wf%Ex, 1), size(wf%Ex, 2), nz]
maxdims = dims
if (grow) maxdims(3) = H5S_UNLIMITED_F

call pmd_init_compound_complex (complex_t)
call H5Pcreate_f (H5P_DATASET_CREATE_F, p_id, h5_err)
if (grow) call H5Pset_chunk_f (p_id, 3, [dims(1), dims(2), 1_hsize_t], h5_err)

call create_component ('x')
if (.not. err .and. allocated(wf%Ey)) call create_component ('y')

call H5Pclose_f (p_id, h5_err)
call pmd_kill_compound_complex (complex_t)
if (err) return

err_flag = .false.

!------------------------------------------------------------------------------
contains

!+
! Subroutine create_component (name)
!
! Routine to create one empty polarisation component dataset with its attributes.
!-

subroutine create_component (name)

character(*) name

!

call H5Screate_simple_f (3, dims, s_id, h5_err, maxdims)
call H5Dcreate_f (m_id, name, complex_t, s_id, d_id, h5_err, p_id)
if (h5_err < 0) then
  call out_io (s_error$, r_name, 'CANNOT CREATE COMPONENT ' // name // ' IN ' // trim(file_name))
  err = .true.
  return
endif
call H5Dclose_f (d_id, h5_err)
call H5Sclose_f (s_id, h5_err)
call write_component_attributes (m_id, name, err)

end subroutine create_component

end subroutine wavefront_openpmd_create

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_openpmd_write_planes (m_id, name, iz, planes, err_flag)
!
! Routine to write a run of planes into one component dataset of a file made by
! wavefront_openpmd_create, at slice index iz onward. An extendable dataset grows to
! hold them. A contiguous one must already hold index iz + size(planes, 3) - 1.
!
! Input:
!   m_id           -- integer(hid_t): The open mesh record group.
!   name           -- character(*): The component, 'x' or 'y'.
!   iz             -- integer: Slice index of the first plane written, from 1.
!   planes(:,:,:)  -- complex(wf_rp): The planes, (nx, ny, n).
!
! Output:
!   err_flag       -- logical: Set True on error, False otherwise.
!-

subroutine wavefront_openpmd_write_planes (m_id, name, iz, planes, err_flag)

complex(wf_rp), target, contiguous :: planes(:,:,:)
character(*) name
integer iz
integer(hid_t) m_id, d_id, complex_t, f_space, m_space
integer(hsize_t) dims(3), start(3), count(3), maxdims(3)
integer h5_err, ndim
logical err_flag
type (c_ptr) f_ptr

!

err_flag = .true.
count = shape(planes)
start = [0_hsize_t, 0_hsize_t, int(iz - 1, hsize_t)]

call H5Dopen_f (m_id, name, d_id, h5_err);  if (h5_err < 0) return
call H5Dget_space_f (d_id, f_space, h5_err)
call H5Sget_simple_extent_dims_f (f_space, dims, maxdims, ndim)
call H5Sclose_f (f_space, h5_err)
if (ndim /= 3 .or. any(dims(1:2) /= count(1:2))) then
  call H5Dclose_f (d_id, h5_err)
  return
endif

if (dims(3) < start(3) + count(3)) then
  dims(3) = start(3) + count(3)
  call H5Dset_extent_f (d_id, dims, h5_err)
  if (h5_err < 0) then
    call H5Dclose_f (d_id, h5_err)
    return
  endif
endif

call H5Dget_space_f (d_id, f_space, h5_err)
call H5Sselect_hyperslab_f (f_space, H5S_SELECT_SET_F, start, count, h5_err)
call H5Screate_simple_f (3, count, m_space, h5_err)
call pmd_init_compound_complex (complex_t)
f_ptr = c_loc(planes)
call H5Dwrite_f (d_id, complex_t, f_ptr, h5_err, m_space, f_space)
call pmd_kill_compound_complex (complex_t)
call H5Sclose_f (m_space, h5_err)
call H5Sclose_f (f_space, h5_err)
call H5Dclose_f (d_id, h5_err)
if (h5_err < 0) return

err_flag = .false.

end subroutine wavefront_openpmd_write_planes

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_openpmd_read_planes (m_id, name, iz, planes, err_flag)
!
! Routine to read a run of planes from one component dataset, from slice index iz
! onward, into planes, whose shape says how many. The inverse of
! wavefront_openpmd_write_planes.
!
! Input:
!   m_id           -- integer(hid_t): The open mesh record group.
!   name           -- character(*): The component, 'x' or 'y'.
!   iz             -- integer: Slice index of the first plane read, from 1.
!
! Output:
!   planes(:,:,:)  -- complex(wf_rp): The planes, (nx, ny, n).
!   err_flag       -- logical: Set True on error, False otherwise.
!-

subroutine wavefront_openpmd_read_planes (m_id, name, iz, planes, err_flag)

complex(wf_rp), target, contiguous :: planes(:,:,:)
character(*) name
integer iz
integer(hid_t) m_id, d_id, complex_t, f_space, m_space
integer(hsize_t) dims(3), start(3), count(3), maxdims(3)
integer h5_err, ndim
logical err_flag
type (c_ptr) f_ptr

!

err_flag = .true.
count = shape(planes)
start = [0_hsize_t, 0_hsize_t, int(iz - 1, hsize_t)]

call H5Dopen_f (m_id, name, d_id, h5_err);  if (h5_err < 0) return
call H5Dget_space_f (d_id, f_space, h5_err)
call H5Sget_simple_extent_dims_f (f_space, dims, maxdims, ndim)
if (ndim /= 3 .or. any(dims(1:2) /= count(1:2)) .or. dims(3) < start(3) + count(3)) then
  call H5Sclose_f (f_space, h5_err)
  call H5Dclose_f (d_id, h5_err)
  return
endif

call H5Sselect_hyperslab_f (f_space, H5S_SELECT_SET_F, start, count, h5_err)
call H5Screate_simple_f (3, count, m_space, h5_err)
call pmd_init_compound_complex (complex_t)
f_ptr = c_loc(planes)
call H5Dread_f (d_id, complex_t, f_ptr, h5_err, m_space, f_space)
call pmd_kill_compound_complex (complex_t)
call H5Sclose_f (m_space, h5_err)
call H5Sclose_f (f_space, h5_err)
call H5Dclose_f (d_id, h5_err)
if (h5_err < 0) return

err_flag = .false.

end subroutine wavefront_openpmd_read_planes

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_openpmd_mesh (f_id, m_id, err_flag)
!
! Routine to open the mesh record group of an open EXT_Wavefront file, the parent of
! the component datasets wavefront_openpmd_write_planes and wavefront_openpmd_read_planes
! take. The caller closes m_id.
!
! Input:
!   f_id      -- integer(hid_t): The open file.
!
! Output:
!   m_id      -- integer(hid_t): The open mesh record group.
!   err_flag  -- logical: Set True when the file has no mesh record, False otherwise.
!-

subroutine wavefront_openpmd_mesh (f_id, m_id, err_flag)

integer(hid_t) f_id, m_id
logical err_flag, err

!

err_flag = .true.
if (.not. hdf5_exists(f_id, trim(mesh_path), err, .false.)) return
m_id = hdf5_open_group (f_id, trim(mesh_path), err, .true.);  if (err) return
err_flag = .false.

end subroutine wavefront_openpmd_mesh

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine write_mesh_attributes (m_id, wf, z_offset, z_coordinate, err_flag)
!
! Routine to write the mesh record's attributes: the base standard's mesh attributes
! and the Wavefront extension's required set (module header), from the grid and the
! wavelength of wf. Shared by the whole-file writer and wavefront_openpmd_create, so the
! two cannot disagree.
!
! Input:
!   m_id          -- integer(hid_t): The open mesh record group.
!   wf            -- wavefront_struct: The grid and the wavelength.
!   z_offset      -- real(rp): The z of gridGlobalOffset [m].
!   z_coordinate  -- real(rp): The beamline position written as zCoordinate [m].
!
! Output:
!   err_flag      -- logical: Set True on error, False otherwise.
!-

subroutine write_mesh_attributes (m_id, wf, z_offset, z_coordinate, err_flag)

type (wavefront_struct) wf
real(rp) z_offset, z_coordinate, gmax_x, gmax_y, e_photon
integer(hid_t) m_id
integer i
logical err_flag, err, bad
! The dimension exponents of one grid axis: a length, so L = 1 and the rest zero.
! gridUnitDimension repeats them once per axis (openPMD 2.0, "Mesh Based Records").
real(rp), parameter :: len_dim(7) = [1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]

!

gmax_x = (size(wf%Ex,1) - 1) * wf%dx / 2
gmax_y = (size(wf%Ex,2) - 1) * wf%dy / 2
e_photon = h_planck * c_light / wf%wavelength * e_charge     ! [J]. h_planck is eV s.

! Each writer reports its own failure, and bad accumulates them, as fel_h5_mod's do.

bad = .false.
call hdf5_write_attribute_string (m_id, 'geometry', 'cartesian', err)
bad = bad .or. err
call hdf5_write_attribute_string (m_id, 'axisLabels', [character(1):: 'z', 'y', 'x'], err)
bad = bad .or. err
call hdf5_write_attribute_real (m_id, 'gridSpacing', [wf%dz, wf%dy, wf%dx], err)
bad = bad .or. err
call hdf5_write_attribute_real (m_id, 'gridGlobalOffset', [z_offset, -gmax_y, -gmax_x], err)
bad = bad .or. err
call hdf5_write_attribute_real (m_id, 'gridUnitSI', [1.0_rp, 1.0_rp, 1.0_rp], err)
bad = bad .or. err
call hdf5_write_attribute_real (m_id, 'gridUnitDimension', [(len_dim, i = 1, 3)], err)
bad = bad .or. err
call hdf5_write_attribute_real (m_id, 'unitDimension', &
                                [1.0_rp, 1.0_rp, -3.0_rp, -1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], err)
bad = bad .or. err
call hdf5_write_attribute_real (m_id, 'timeOffset', 0.0_rp, err)
bad = bad .or. err
call hdf5_write_attribute_real (m_id, 'photonEnergy', e_photon, err)
bad = bad .or. err
call hdf5_write_attribute_string (m_id, 'temporalDomain', 'time', err)
bad = bad .or. err
call hdf5_write_attribute_string (m_id, 'spatialDomain', 'r', err)
bad = bad .or. err
call hdf5_write_attribute_real (m_id, 'zCoordinate', z_coordinate, err)
err_flag = bad .or. err

end subroutine write_mesh_attributes

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine write_component_attributes (m_id, name, err_flag)
!
! Routine to write the attributes of one polarisation component dataset: the units,
! unitSI 1 since the field is in V/m, and the required openPMD attribute "position",
! the sample's offset within its cell in units of gridSpacing. A wavefront's samples
! sit on the cell corners, so it is zero on all three axes.
!
! Note: The attributes go on the dataset. So this routine takes the dataset name and
! hands it to the H5LT calls that write through a parent. The plain attribute writers
! take an object id. Passing them a path as the attribute name puts an attribute with
! slashes in its name on the root group instead.
!
! Input:
!   m_id      -- integer(hid_t): The open mesh record group, the dataset's parent.
!   name      -- character(*): The component, 'x' or 'y'.
!
! Output:
!   err_flag  -- logical: Set True on error, False otherwise.
!-

subroutine write_component_attributes (m_id, name, err_flag)

character(*) name
integer(hid_t) m_id
integer h5_err
logical err_flag, err

!

call pmd_write_units_to_dataset (m_id, name, name, unit_V_per_m, err)
call H5LTset_attribute_double_f (m_id, name, 'position', [0.0_rp, 0.0_rp, 0.0_rp], &
                                 3_size_t, h5_err)
err_flag = err .or. (h5_err < 0)

end subroutine write_component_attributes

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_read_openpmd (wf, file_name, err_flag, photon_energy)
!
! Routine to read an openPMD EXT_Wavefront file in the layout of the module header.
! Anything this reader cannot represent is refused: a frequency-domain or k-space
! field, an axis order other than (z,y,x), a missing required attribute.
!
! photon_energy is returned so the caller can match the file to the field set entry
! whose harmonic it carries. wf%wavelength is derived from it.
!
! Input:
!   file_name      -- character(*): File to read.
!
! Output:
!   wf             -- wavefront_struct: Wavefront read from the file.
!   err_flag       -- logical: Set True on error, False otherwise.
!   photon_energy  -- real(rp): Central photon energy from the file [J].
!-

subroutine wavefront_read_openpmd (wf, file_name, err_flag, photon_energy)

type (wavefront_struct), target :: wf
character(*) file_name
real(rp) photon_energy
real(rp) spacing(3), offset(3), gunit(3), rval
integer(hid_t) f_id, m_id, complex_t
integer h5_err, ndim(3), ndim_y(3)
logical err_flag, err
character(40) sval
character(4), allocatable :: labels(:)
type (hdf5_info_struct) info
character(*), parameter :: r_name = 'wavefront_read_openpmd'

!

err_flag = .true.

call hdf5_open_file (file_name, 'READ', f_id, err);  if (err) return

if (.not. hdf5_exists(f_id, trim(mesh_path), err, .false.)) then
  call out_io (s_error$, r_name, 'NO ' // trim(mesh_path) // ' MESH RECORD IN ' // trim(file_name), &
                                 '(EXT_Wavefront names the record electricField).')
  return
endif
m_id = hdf5_open_group (f_id, trim(mesh_path), err, .true.);  if (err) return

! The Wavefront extension's required attributes. Each missing one is named.

call hdf5_read_attribute_string (m_id, 'temporalDomain', sval, err, .false.)
if (err) then
  call out_io (s_error$, r_name, 'REQUIRED ATTRIBUTE temporalDomain MISSING IN ' // trim(file_name))
  return
endif
if (sval /= 'time') then
  call out_io (s_error$, r_name, 'temporalDomain "' // trim(sval) // '" IN ' // trim(file_name), &
                                 'ONLY THE time DOMAIN (field in V/m) IS IMPLEMENTED.')
  return
endif

call hdf5_read_attribute_string (m_id, 'spatialDomain', sval, err, .false.)
if (err) then
  call out_io (s_error$, r_name, 'REQUIRED ATTRIBUTE spatialDomain MISSING IN ' // trim(file_name))
  return
endif
if (sval /= 'r') then
  call out_io (s_error$, r_name, 'spatialDomain "' // trim(sval) // '" IN ' // trim(file_name), &
                                 'ONLY CARTESIAN r SPACE IS IMPLEMENTED.')
  return
endif

call hdf5_read_attribute_real (m_id, 'photonEnergy', photon_energy, err, .false.)
if (err) then
  call out_io (s_error$, r_name, 'REQUIRED ATTRIBUTE photonEnergy MISSING IN ' // trim(file_name))
  return
endif

call hdf5_read_attribute_real (m_id, 'zCoordinate', rval, err, .false.)
if (err) then
  call out_io (s_error$, r_name, 'REQUIRED ATTRIBUTE zCoordinate MISSING IN ' // trim(file_name))
  return
endif
wf%s_position = rval

call hdf5_read_attribute_string (m_id, 'axisLabels', labels, err, .false.)
if (err) then
  call out_io (s_error$, r_name, 'REQUIRED ATTRIBUTE axisLabels MISSING IN ' // trim(file_name))
  return
endif
if (size(labels) /= 3) then
  call out_io (s_error$, r_name, 'axisLabels MUST HAVE THREE ENTRIES IN ' // trim(file_name))
  return
endif
if (labels(1) /= 'z' .or. labels(2) /= 'y' .or. labels(3) /= 'x') then
  call out_io (s_error$, r_name, 'axisLabels (' // trim(labels(1)) // ',' // trim(labels(2)) // &
               ',' // trim(labels(3)) // ') IN ' // trim(file_name), &
               'ONLY THE (z, y, x) STORED ORDER IS IMPLEMENTED.')
  return
endif

call hdf5_read_attribute_real (m_id, 'gridSpacing', spacing, err, .false.)
if (err) then
  call out_io (s_error$, r_name, 'REQUIRED ATTRIBUTE gridSpacing MISSING IN ' // trim(file_name))
  return
endif
! openPMD states the grid in its own units with the SI factor beside them. A file
! written elsewhere may use units other than metres, so the spacing is scaled by it.

call hdf5_read_attribute_real (m_id, 'gridUnitSI', gunit, err, .false.)
if (err) then
  call out_io (s_error$, r_name, 'REQUIRED ATTRIBUTE gridUnitSI MISSING IN ' // trim(file_name))
  return
endif
if (any(gunit <= 0)) then
  call out_io (s_error$, r_name, 'gridUnitSI MUST BE POSITIVE ON EVERY AXIS IN ' // trim(file_name))
  return
endif
spacing = spacing * gunit
wf%dz = spacing(1);  wf%dy = spacing(2);  wf%dx = spacing(3)
wf%wavelength = h_planck * c_light * e_charge / photon_energy

call hdf5_read_attribute_real (m_id, 'gridGlobalOffset', offset, err, .false.)   ! Read but unused.
                                            ! The walk's grid is centered by convention.

! The components. Dataset dims come back in Fortran order (nx, ny, nz), by the same API
! symmetry the writer used. Both components are read into buffers of the shape of x.
! An HDF5 read has no bound of its own, so y must be shaped as x is. A larger y would
! write past the end of the buffer.

info = hdf5_object_info (m_id, 'x', err, .true.);  if (err) return
ndim = int(info%data_dim(1:3))
if (any(ndim <= 0)) then
  call out_io (s_error$, r_name, 'COMPONENT x IS NOT A THREE-DIMENSIONAL DATASET IN ' // trim(file_name))
  return
endif
if (allocated(wf%Ex)) deallocate (wf%Ex)
if (allocated(wf%Ey)) deallocate (wf%Ey)
allocate (wf%Ex(ndim(1), ndim(2), ndim(3)))

call pmd_init_compound_complex (complex_t)
call read_component ('x', wf%Ex, err)
if (.not. err .and. hdf5_exists(m_id, 'y', err, .false.)) then
  info = hdf5_object_info (m_id, 'y', err, .true.)
  if (.not. err) then
    ndim_y = int(info%data_dim(1:3))
    if (any(ndim_y /= ndim)) then
      call out_io (s_error$, r_name, 'COMPONENT y IS SHAPED (\i0\, \i0\, \i0\) WHERE x IS ' // &
                   '(\i0\, \i0\, \i0\) IN ' // trim(file_name), &
                   'BOTH POLARIZATIONS OF ONE FIELD SHARE ONE GRID.', &
                   i_array = [ndim_y(3), ndim_y(2), ndim_y(1), ndim(3), ndim(2), ndim(1)])
      err = .true.
    endif
  endif
  if (.not. err) then
    allocate (wf%Ey(ndim(1), ndim(2), ndim(3)))
    call read_component ('y', wf%Ey, err)
  endif
endif
call pmd_kill_compound_complex (complex_t)
if (err) return

call h5gclose_f (m_id, h5_err)
call h5fclose_f (f_id, h5_err)

err_flag = .false.

!------------------------------------------------------------------------------
contains

!+
! Subroutine read_component (name, fld, cerr)
!
! Routine to read one polarisation component dataset of the mesh record into fld.
!-

subroutine read_component (name, fld, cerr)

character(*) name
complex(wf_rp), target, contiguous :: fld(:,:,:)
real(rp) unit
logical cerr, uerr
integer(hid_t) d_id
type(c_ptr) f_ptr

!

cerr = .true.
call h5dopen_f (m_id, name, d_id, h5_err);  if (h5_err < 0) return

! The samples are in the file's own unit and unitSI takes them to V/m. openPMD requires
! the attribute on every component, and this writer puts 1 there.

call hdf5_read_attribute_real (d_id, 'unitSI', unit, uerr, .false.)
if (uerr) then
  call out_io (s_error$, r_name, 'REQUIRED ATTRIBUTE unitSI MISSING ON COMPONENT ' // trim(name) // &
                                 ' IN ' // trim(file_name))
  call h5dclose_f (d_id, h5_err);  return
endif
if (unit <= 0) then
  call out_io (s_error$, r_name, 'unitSI MUST BE POSITIVE ON COMPONENT ' // trim(name) // ' IN ' // trim(file_name))
  call h5dclose_f (d_id, h5_err);  return
endif

f_ptr = c_loc(fld)
call h5dread_f (d_id, complex_t, f_ptr, h5_err)
call h5dclose_f (d_id, h5_err)
if (h5_err < 0) return
if (unit /= 1) fld = fld * real(unit, wf_rp)
cerr = .false.

end subroutine read_component

end subroutine wavefront_read_openpmd

end module wavefront_openpmd_mod
