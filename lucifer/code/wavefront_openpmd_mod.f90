!+
! Module wavefront_openpmd_mod
!
! openPMD EXT_Wavefront I/O for wavefront_struct: the openPMD-standard base plus the
! Wavefront extension (openPMD-standard, branch upcoming-2.0.0, EXT_Wavefront.md).
! The standard document is authoritative for this format. This module and the
! harness's h5py reader validate against its text independently.
!
! Layout decisions (see the physics manual, fel-physics.md sec-field-set):
!
!   /                        openPMD "2.0.0", openPMDextension "Wavefront",
!                            basePath "/data/%T/", meshesPath "meshes/",
!                            iterationEncoding "groupBased", iterationFormat "/data/%T/"
!   /data/1/                 one iteration per file; time = 0, dt = 0, timeUnitSI = 1
!   /data/1/meshes/electricField
!                            the one mesh record. Attributes on the mesh record (the
!                            extension's heading. Its body says "series", and the
!                            contradiction is resolved here in the record's favor:
!                            photonEnergy is a property of one field):
!                            geometry 'cartesian', axisLabels ['z','y','x'],
!                            gridSpacing [dz,dy,dx], gridGlobalOffset [0,-gmax,-gmax],
!                            gridUnitSI [1,1,1], gridUnitDimension (a length per axis),
!                            unitDimension (1,1,-3,-1,0,0,0) (V/m),
!                            timeOffset 0, photonEnergy [J], temporalDomain 'time',
!                            spatialDomain 'r', zCoordinate [m].
!   .../electricField/x      complex compound {r,i} dataset, which h5py reads natively.
!                            Component attributes: unitSI 1 and position [0,0,0], the
!                            sample's offset within its cell in units of gridSpacing.
!   .../electricField/y      present only when the wavefront carries Ey. Both transverse
!                            polarizations live in one file as components, the
!                            improvement over the Genesis format's one-per-file. The z
!                            component is never written: a paraxial code has none, and
!                            an absent component is ordinary openPMD.
!
! The dataset is stored exactly as the Fortran (nx, ny, nslice) array, which the HDF5
! Fortran API records as a C-order (nslice, ny, nx) dataspace: zero-copy, and
! numpy-natural for per-slice access. axisLabels declare that stored order. The slice
! axis is the one labeled z (third in the logical x,y,z reading). Slices are
! simultaneous, so they are a mesh axis, never the openPMD iteration.
!
! Reading implements exactly what writing produces and refuses the rest:
! temporalDomain 'frequency', spatialDomain 'k', an axisLabels order other than
! (z,y,x), and any missing required attribute.
!-

module wavefront_openpmd_mod

use wavefront_mod
use hdf5_openpmd_mod
use sim_utils

implicit none

character(*), parameter, private :: iter_path = 'data/1'
character(*), parameter, private :: mesh_path = 'data/1/meshes/electricField'

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function file_is_openpmd (file_name) result (is_pmd)
!
! Routine to test whether a file is an openPMD file. This is what the import path refuses
! on, and it holds for any openPMD file, particles as well as fields: an openPMD file
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
! Subroutine wavefront_write_openpmd (wf, file_name, s_pos, err_flag)
!
! Routine to write one wavefront as an openPMD EXT_Wavefront file, using the layout
! described in the module header. The field is in V/m already (temporalDomain 'time'),
! so unitSI = 1 throughout.
!
! Input:
!   wf         -- wavefront_struct: The field. Ey written when allocated.
!   file_name  -- character(*): Output file.
!   s_pos      -- real(rp): Lattice position of the dump plane [m] (zCoordinate).
!
! Output:
!   err_flag   -- logical: Set True on error, False otherwise.
!-

subroutine wavefront_write_openpmd (wf, file_name, s_pos, err_flag)

type (wavefront_struct), target :: wf
character(*) file_name
real(rp) s_pos, gmax_x, gmax_y, e_photon
integer(hid_t) f_id, it_id, m_id, complex_t
integer h5_err, i
logical err_flag, err
character(*), parameter :: r_name = 'wavefront_write_openpmd'
! The dimension exponents of one grid axis: a length, so L = 1 and the rest zero.
! gridUnitDimension repeats them once per axis (openPMD 2.0, "Mesh Based Records").
real(rp), parameter :: len_dim(7) = [1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]

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

gmax_x = (size(wf%Ex,1) - 1) * wf%dx / 2
gmax_y = (size(wf%Ex,2) - 1) * wf%dy / 2
e_photon = h_planck * c_light / wf%wavelength * e_charge     ! [J]. h_planck is eV s.

call hdf5_write_attribute_string (m_id, 'geometry', 'cartesian', err)
call hdf5_write_attribute_string (m_id, 'axisLabels', [character(1):: 'z', 'y', 'x'], err)
call hdf5_write_attribute_real (m_id, 'gridSpacing', [wf%dz, wf%dy, wf%dx], err)
call hdf5_write_attribute_real (m_id, 'gridGlobalOffset', [0.0_rp, -gmax_y, -gmax_x], err)
call hdf5_write_attribute_real (m_id, 'gridUnitSI', [1.0_rp, 1.0_rp, 1.0_rp], err)
call hdf5_write_attribute_real (m_id, 'gridUnitDimension', [(len_dim, i = 1, 3)], err)
call hdf5_write_attribute_real (m_id, 'unitDimension', [1.0_rp, 1.0_rp, -3.0_rp, -1.0_rp, 0.0_rp, 0.0_rp, 0.0_rp], err)
call hdf5_write_attribute_real (m_id, 'timeOffset', 0.0_rp, err)
call hdf5_write_attribute_real (m_id, 'photonEnergy', e_photon, err)
call hdf5_write_attribute_string (m_id, 'temporalDomain', 'time', err)
call hdf5_write_attribute_string (m_id, 'spatialDomain', 'r', err)
call hdf5_write_attribute_real (m_id, 'zCoordinate', s_pos, err)

call pmd_init_compound_complex (complex_t)

call pmd_write_complex_to_dataset (m_id, 'x', complex_t, 'x', unit_V_per_m, wf%Ex, err)
if (err) then
  call pmd_kill_compound_complex (complex_t)
  return
endif
call write_position ('x')

if (allocated(wf%Ey)) then
  call pmd_write_complex_to_dataset (m_id, 'y', complex_t, 'y', unit_V_per_m, wf%Ey, err)
  if (err) then
    call pmd_kill_compound_complex (complex_t)
    return
  endif
  call write_position ('y')
endif

call pmd_kill_compound_complex (complex_t)
call h5gclose_f (m_id, h5_err)
call h5fclose_f (f_id, h5_err)

err_flag = .false.

!------------------------------------------------------------------------------
contains

!+
! Subroutine write_position (name)
!
! Routine to write the required openPMD component attribute "position" on one
! polarisation component dataset. The value is the sample's offset within its cell in
! units of gridSpacing, and a wavefront's samples sit on the cell corners, so it is zero
! on all three axes.
!
! The attribute goes on the dataset, which is why this takes the dataset name and hands it
! to the H5LT call that writes through a parent: the plain attribute writers take an object
! id, and passing a path as the attribute name puts a slash-laden attribute on the root
! group instead. Caught by the validation's file comparison against the Python writer.
!-

subroutine write_position (name)

character(*) name

!

call H5LTset_attribute_double_f (m_id, name, 'position', [0.0_rp, 0.0_rp, 0.0_rp], 3_size_t, h5_err)

end subroutine write_position

end subroutine wavefront_write_openpmd

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_read_openpmd (wf, file_name, err_flag, photon_energy)
!
! Routine to read an openPMD EXT_Wavefront file written to the module header's layout.
! Everything this reader cannot represent is refused: a frequency-domain or
! k-space field, an axis order other than (z,y,x), a missing required attribute.
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
real(rp) spacing(3), offset(3), rval
integer(hid_t) f_id, m_id, complex_t
integer h5_err, ndim(3)
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
wf%ref_position = rval

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
wf%dz = spacing(1);  wf%dy = spacing(2);  wf%dx = spacing(3)
wf%wavelength = h_planck * c_light * e_charge / photon_energy

call hdf5_read_attribute_real (m_id, 'gridGlobalOffset', offset, err, .false.)   ! Read, unused:
                                            ! the walk's grid is centered by convention.

! The components. Dataset dims come back in Fortran order (nx, ny, nz), the same API
! symmetry the writer used.

info = hdf5_object_info (m_id, 'x', err, .true.);  if (err) return
ndim = int(info%data_dim(1:3))
if (allocated(wf%Ex)) deallocate (wf%Ex)
if (allocated(wf%Ey)) deallocate (wf%Ey)
allocate (wf%Ex(ndim(1), ndim(2), ndim(3)))

call pmd_init_compound_complex (complex_t)
call read_component ('x', wf%Ex, err)
if (.not. err .and. hdf5_exists(m_id, 'y', err, .false.)) then
  allocate (wf%Ey(ndim(1), ndim(2), ndim(3)))
  call read_component ('y', wf%Ey, err)
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
logical cerr
integer(hid_t) d_id
type(c_ptr) f_ptr

!

cerr = .true.
call h5dopen_f (m_id, name, d_id, h5_err);  if (h5_err < 0) return
f_ptr = c_loc(fld)
call h5dread_f (d_id, complex_t, f_ptr, h5_err)
call h5dclose_f (d_id, h5_err)
if (h5_err < 0) return
cerr = .false.

end subroutine read_component

end subroutine wavefront_read_openpmd

end module wavefront_openpmd_mod
