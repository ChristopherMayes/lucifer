!+
! Module fel_io_mod
!
! Routines to write the file artifacts of a run: the field-set dumps (Genesis and openPMD
! formats), the openPMD field imports, the escaped-slice bank and its drain, the
! full-pulse reconstruction, the stats-file finalize, and the wake-eloss block writer.
!
! Every routine takes the run state (fel_run_struct) explicitly. Nothing here stops.
! Errors return through err_flag and the caller decides (docs/user-guide.md). All
! terminal output goes through out_io. The file-name echo lines are s_blank$, so
! scripts see bare names.
!-

module fel_io_mod

use fel_struct
use fel_input_mod
use fel_schema_mod
use fel_timer_mod
use fel_cost_mod
use fel_print_mod
use fel_record_mod
use wavefront_openpmd_mod
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
! For bp_com%num_lat_files, the parser's tally of how many files it opened. Read only in
! fel_write_meta.
use bmad_parser_struct, only: bp_com
! For bmad_version_date, the distribution's own version. Tao reports the same string.
use bmad_version_mod, only: bmad_version_date
! For lucifer_version$, the version of this build, written by CMake into the build directory.
use fel_version_mod, only: lucifer_version$
!$ use omp_lib

implicit none

! The keys of the noise, which the header writes in a part of their own, the defaults it
! shows beside the keys the input set, and the units that take a prefix.

character(32), parameter, private :: noise_keys(4) = [character(32) :: &
      'quiet_start', 'shot_noise', 'shot_noise_method', 'beamlet_min_electrons']
character(48), parameter, private :: shown_defaults(14) = [character(48) :: &
      'beam.load_mode', 'beam.quiet_start', 'beam.shot_noise', 'beam.shot_noise_method', &
      'field.lambda0', 'field.harmonics', 'field.seed_power', 'run.out_root', 'run.ran_seed', &
      'run.source_filter', 'run.migrate', 'bmad.sr_wakes_on', 'bmad.csr_and_space_charge_on', &
      'bmad.radiation_damping_on']
character(16), parameter, private :: scaled_units(9) = [character(16) :: &
      'm', 'C', 'W', 'eV', 'rad', 'A', 'J', 's', 'm rad']

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_dump_beam (run, ele, prefix, err_flag, s_pos)
!
! Routine to write the beam at the given filename prefix as openPMD (.beam.h5). This is
! the only particle format this tracker writes. Genesis .par.h5 conversion lives in
! lucifer/tests/scripts/convert_genesis.py, outside the physics.
!
! Input:
!   run       -- fel_run_struct: Run state.
!   ele       -- ele_struct: Element the beam sits at.
!   prefix    -- character(*): Filename prefix. Format suffixes are appended.
!   s_pos     -- real(rp), optional: Where along ele the beam sits. A dump at an element's
!                  end passes ele%s. The particle records then carry the downstream face,
!                  and ele is the element completed. If not present, the beam sits at the
!                  upstream face (the initial dump), and the element completed is the one
!                  before ele.
!
! Output:
!   err_flag  -- logical: Set True if a file could not be written. False otherwise.
!-

subroutine fel_dump_beam (run, ele, prefix, err_flag, s_pos)

type (fel_run_struct), target :: run
type (ele_struct) ele
type (ele_struct), pointer :: done
real(rp), optional :: s_pos
integer ix_done
character(*) prefix
logical err_flag
logical eerr

!

err_flag = .true.

call fel_write_openpmd_beam (run%fbeam, ele, prefix // '.beam.h5', eerr, s_pos = s_pos)
if (eerr) return

! The checkpoint group, where this boundary is one a continuation may start from.

ix_done = ele%ix_ele
if (.not. present(s_pos)) ix_done = ele%ix_ele - 1
if (fel_checkpoint_eligible(run, ix_done)) then
  done => run%lat%branch(0)%ele(ix_done)
  call fel_write_checkpoint_group (prefix // '.beam.h5', done%s, ix_done, done%name, eerr, run%fbeam)
  if (eerr) return
endif

err_flag = .false.

end subroutine fel_dump_beam

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_checkpoint_eligible (run, ix_done) result (ok)
!
! Routine to say whether the boundary at the end of element ix_done is one a continuation
! may start from. Frames written there then carry the checkpoint group
! (fel_write_checkpoint_group). The boundary is eligible outside every undulator and at
! the physical end of one, with the beam in the averaged chart.
!
! A sliced endpoint inside an undulator is not a completed handoff. There a superimposed
! element has cut one physical device into super_slaves. The averaged tracker's px is
! the slow momentum inside a field that has not ended. The unaveraged tracker's ramp has
! closed on a plane the device does not end at. The frame is written whatever the
! answer. Only the group is withheld, and a continuation from that frame is refused.
!
! Note: A lattice with a cut undulator is refused at setup, since its slaves refer to
! their lord for the field model (fel_assert_wiggler_sane). So the slave branch below
! guards a boundary no run reaches yet.
!
! Input:
!   run      -- fel_run_struct: Run state.
!   ix_done  -- integer: Index of the element completed at the boundary, 0 for the entry face.
!
! Output:
!   ok       -- logical: True if a checkpoint is written at this boundary.
!-

function fel_checkpoint_eligible (run, ix_done) result (ok)

type (fel_run_struct), target :: run
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele, lord
real(rp) s_b, s1, s2, tol
integer ix_done, je, il
logical ok

!

ok = .false.
if (run%fbeam%quiver_in_px) return
branch => run%lat%branch(0)
if (ix_done < 0 .or. ix_done > branch%n_ele_track) return
s_b = branch%ele(ix_done)%s

! Refuse a boundary strictly inside an undulator's physical extent. The extent of a
! slave is its lord's, since the lord is the device. The tolerance (1 nm) is a length
! no lattice resolves.

tol = 1.0e-9_rp
do je = 1, branch%n_ele_track
  if (.not. run%is_fel(je)) cycle
  ele => branch%ele(je)
  s1 = ele%s_start
  s2 = ele%s
  if (ele%slave_status == super_slave$) then
    do il = 1, ele%n_lord
      lord => pointer_to_lord(ele, il)
      if (lord%key /= ele%key) cycle
      s1 = lord%s_start
      s2 = lord%s
      exit
    enddo
  endif
  if (s_b > s1 + tol .and. s_b < s2 - tol) return
enddo

ok = .true.

end function fel_checkpoint_eligible

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_frame_attributes (file_name, run, ie, err_flag)
!
! Routine to stamp a dump with where it was taken.
!
! A frame is read on its own, so it has to say where along the line it is and what the
! beam was moving through. The particle records are the coordinates the selected
! tracking method produced there (docs/reading-output.md). The period-averaged map has
! no fast oscillation, and the unaveraged mode has it resolved. So the method and the
! device's numbers go in the file, and the reader needs no lattice. A frame taken in a
! break carries the element and no undulator numbers.
!
! The attributes go under lucifer/frame (fel_private_group) in a file the writers have
! already closed. So no writer's layout changes, every kind of file is stamped the same
! way, and the root keeps the standard's declarations. The iteration's time is the
! standard's own and goes where the standard puts it.
!
! Input:
!   file_name -- character(*): An openPMD file the writers have finished.
!   run       -- fel_run_struct: Run state, read for the element and the undulator.
!   ie        -- integer: Index of the element the frame was taken in.
!
! Output:
!   err_flag  -- logical: Set True if the file could not be stamped. False otherwise.
!-

subroutine fel_frame_attributes (file_name, run, ie, err_flag)

type (fel_run_struct), target :: run
type (fel_und_struct), pointer :: und
type (branch_struct), pointer :: branch
type (hdf5_info_struct) obj_info
type (floor_position_struct) fl
integer(hid_t) f_id, g_id, d_id, it_id
integer ie, h5_err
real(rp) t_ref, ds_ele
logical err_flag, err
character(*) file_name
character(200) it_name
character(*), parameter :: r_name = 'fel_frame_attributes'

!

err_flag = .true.
branch => run%lat%branch(0)

call hdf5_open_file (file_name, 'APPEND', f_id, err);  if (err) return
g_id = fel_private_group (f_id, 'frame', .true., err)
if (err) then
  call h5fclose_f (f_id, h5_err)
  return
endif

! Frame format 3.0: the attributes sit under lucifer/frame. A .beam.h5 at this version
! holds the coordinates the selected tracking method produced, wherever it was taken.
call hdf5_write_attribute_string (g_id, 'frameFormat', fel_frame_format$, err)
call hdf5_write_attribute_real (g_id, 'sPosition', run%z_now, err)
call hdf5_write_attribute_string (g_id, 'elementName', trim(branch%ele(ie)%name), err)
call hdf5_write_attribute_int (g_id, 'elementIndex', ie, err)
call hdf5_write_attribute_real (g_id, 'phi0', run%fbeam%phi0, err)
call hdf5_write_attribute_int (g_id, 'sliceFirst', run%dump_is1, err)
call hdf5_write_attribute_int (g_id, 'sliceLast', run%dump_is2, err)

! Where the frame sits in the lab. s is a curvilinear coordinate, so a scene holding a
! line with a bend needs the floor to place the frame and the grid transverse to it.
! The floor is taken at the frame's own s and not the element's, since a frame lands
! mid-element and ele%floor is the downstream face. coords_curvilinear_to_floor walks
! the geometry, so a bend is not straightened.

fl = coords_curvilinear_to_floor ([0.0_rp, 0.0_rp, run%z_now], branch, err)
if (.not. err) then
  call hdf5_write_attribute_real (g_id, 'floorPosition', [fl%r(1), fl%r(2), fl%r(3)], err)
  call hdf5_write_attribute_real (g_id, 'floorAngles', [fl%theta, fl%phi, fl%psi], err)
endif

! openPMD orders a series by the iteration's time and has no notion of s. So the frame
! carries Bmad's reference time rather than s/c. The reference momentum is constant
! through an undulator, so the time inside an element is linear between the element's
! two faces. The writers leave the attribute at zero, which would stack a whole series
! at one instant.

t_ref = branch%ele(ie)%ref_time
if (ie >= 1) then
  ds_ele = branch%ele(ie)%s - branch%ele(ie-1)%s
  if (ds_ele > 0) t_ref = branch%ele(ie-1)%ref_time + (branch%ele(ie)%ref_time - &
                          branch%ele(ie-1)%ref_time) * (run%z_now - branch%ele(ie-1)%s) / ds_ele
endif

! The iteration group's name differs between the two writers, so it is found rather than
! spelled. Each file holds exactly one.

d_id = -1
if (hdf5_exists(f_id, 'data', err, .false.)) d_id = hdf5_open_group (f_id, 'data', err, .true.)
if (d_id >= 0 .and. .not. err) then
  if (hdf5_group_n_links(d_id, err) == 1) then
    call hdf5_get_object_by_index (d_id, 0, it_name, obj_info, err)
    if (.not. err) then
      it_id = hdf5_open_group (d_id, trim(it_name), err, .true.)
      if (.not. err) then
        call hdf5_write_attribute_real (it_id, 'time', t_ref, err)
        call hdf5_write_attribute_real (it_id, 'timeUnitSI', 1.0_rp, err)
        call H5Gclose_f (it_id, h5_err)
      endif
    endif
  endif
  call H5Gclose_f (d_id, h5_err)
endif

! The lattice's elements run from 0, the entry face, and the run's per-element arrays
! run from 1. So the frame taken before the walk enters its first element has a name and
! a position, and no method and no undulator.

if (ie >= 1) then
  call hdf5_write_attribute_string (g_id, 'felMethod', &
        trim(fel_method_name(max(1, min(3, run%fel_mode(ie))))), err)

  ! The undulator, where there is one. A break has an element and no aw.

  if (run%is_fel(ie)) then
    und => run%und_of(ie)
    call hdf5_write_attribute_real (g_id, 'aw', und%aw, err)
    call hdf5_write_attribute_real (g_id, 'ku', und%ku, err)
    call hdf5_write_attribute_real (g_id, 'tilt', und%tilt, err)
    call hdf5_write_attribute_int (g_id, 'helical', merge(1, 0, und%helical), err)

    ! The rest of what rebuilding the orbit takes: where the frame sits in the device, and
    ! how long the device and its ramped ends are. A reader that walks a frame series
    ! gets the lattice's geometry from the frames alone (docs/reading-output.md).

    call hdf5_write_attribute_real (g_id, 'sElement', run%z_now - branch%ele(ie)%s_start, err)
    call hdf5_write_attribute_real (g_id, 'elementLength', branch%ele(ie)%value(l$), err)
    call hdf5_write_attribute_real (g_id, 'rampPeriods', run%fel_ramp(ie), err)
  endif
endif

call H5Gclose_f (g_id, h5_err)
call h5fclose_f (f_id, h5_err)
err_flag = .false.

end subroutine fel_frame_attributes

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_write_slippage_residual (file_name, accuslip, err_flag)
!
! Routine to stamp a field file with the slippage the record has accumulated and not yet
! rotated (fel_slip_struct%accuslip). It is written as the attribute slippageResidual of
! the group lucifer/field (fel_private_group), in fundamental wavelengths, signed.
!
! This is the one piece of the field's state the record itself does not hold. The
! rotation index is folded into the time order the file is written in. This remainder
! decides when the next rotation falls. A restart that starts it at zero rotates on a
! different schedule from the run it continues. The value is written as it stands and
! is never reduced into [0, n_wavelength). The threshold is 0.8 n_wavelength, and a
! residual after a rotation is negative as often as not.
!
! Input:
!   file_name -- character(*): The field file, already written.
!   accuslip  -- real(rp): The record's accumulated slippage [fundamental wavelengths].
!
! Output:
!   err_flag  -- logical: Set True if the attribute could not be written. False otherwise.
!-

subroutine fel_write_slippage_residual (file_name, accuslip, err_flag)

integer(hid_t) f_id, g_id
integer h5_err
real(rp) accuslip
logical err_flag, err
character(*) file_name

!

err_flag = .true.
call hdf5_open_file (file_name, 'APPEND', f_id, err);  if (err) return
g_id = fel_private_group (f_id, 'field', .true., err)
if (err) then
  call h5fclose_f (f_id, h5_err)
  return
endif
call hdf5_write_attribute_real (g_id, 'slippageResidual', accuslip, err)
call H5Gclose_f (g_id, h5_err)
call h5fclose_f (f_id, h5_err)
err_flag = err

end subroutine fel_write_slippage_residual

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_dump_frame (run, ie, at_end, err_flag)
!
! Routine to write one frame of the series run.dump_at_comb asks for: the beam and
! the field set at this comb position. It uses the same writers as the element-end
! dumps. The files are <out_root>-<record>.beam.h5 and <out_root>-<record>.wf.h5, where
! record is the stats row this frame sits on. So a frame and its row share one index.
! The beam frame holds the coordinates the selected tracking method produced there,
! inside an averaged undulator as anywhere else, with nothing reconstructed.
!
! The field's records are written in time order. The whole window of an averaged element
! is rotated in place, written and rotated back, so the run continues from the state it
! had. fel_dump_field_set leaves the records unrotated. That is right at the end of a run
! and wrong in the middle of one. A diagnostic that changes the state it observes is the
! failure the FP32 twin's read-only proof exists to catch. Every other field frame is
! written from a scratch copy of the slices it writes and no others. Each slice is taken
! from its record through fel_field_index. Inside an unaveraged segment the copy is
! carried back half a substep onto the frame's plane, one transform pair a slice with the
! slices in parallel. So a slice-range frame costs the slices it writes. The live record
! is not touched on either path.
!
! Input:
!   run       -- fel_run_struct: Run state.
!   ie        -- integer: Index of the element the frame is taken in.
!   at_end    -- logical: The frame is at the element's end. The walk's z_now reaches an
!                  end by accumulation and can miss ele%s in its last bit. So this flag,
!                  and not a comparison of the two, places the beam at the downstream
!                  face. A frame at an element end is then exactly the element's own dump.
!
! Output:
!   err_flag  -- logical: Set True if a file could not be written. False otherwise.
!-

subroutine fel_dump_frame (run, ie, at_end, err_flag)

type (fel_run_struct), target :: run
type (fel_field_struct), pointer :: ffield(:)
type (fel_beam_struct), pointer :: fbeam
type (wavefront_struct), pointer :: wf
type (wavefront_struct) wf_sub
real(rp) s_frame
integer ie, ihh, first_was, nsub, nslice, j, ifld, js1, js2, m
logical at_end, err_flag, eerr, any_err, whole, ckpt
character(200) prefix
character(8) hsuf

!

err_flag = .true.
ffield => run%ffield

write (prefix, '(2a, i6.6)') trim(run%global%out_root), '-', run%stats%irec

! A frame is written wherever the comb falls, as the beam stands there. Inside an
! unaveraged segment px is the resolved kinetic momentum. Inside an averaged one it is
! the map's period-averaged momentum. The frame carries whichever the selected method
! produced, with nothing reconstructed. fel_dump_beam, whose dumps land on element
! boundaries, keeps the chart check. A frame at an element end that carries the whole
! window is a checkpoint where the boundary is eligible (fel_checkpoint_eligible). It
! gets the group on both of its files. Every other frame is written without it, and a
! continuation from it is refused.

s_frame = run%z_now
if (at_end .and. ie >= 1) s_frame = run%lat%branch(0)%ele(ie)%s
whole = (run%dump_is1 == 1 .and. run%dump_is2 == run%nslice)
ckpt = at_end .and. whole
if (ckpt) ckpt = fel_checkpoint_eligible(run, ie)

! The walk's own beam is written as it stands and never touched. A diagnostic that moved
! the state it observes would steer the run. A run with frames on is the run with them
! off, to the bit.

fbeam => run%fbeam
call fel_write_openpmd_beam (fbeam, run%lat%branch(0)%ele(ie), &
                             trim(prefix) // '.beam.h5', eerr, run%dump_is1, run%dump_is2, &
                             kinetic_ok = .true., s_pos = s_frame)
if (eerr) return
call fel_frame_attributes (trim(prefix) // '.beam.h5', run, ie, eerr)
if (eerr) return
if (ckpt) then
  call fel_write_checkpoint_group (trim(prefix) // '.beam.h5', run%lat%branch(0)%ele(ie)%s, ie, &
                                   run%lat%branch(0)%ele(ie)%name, eerr, run%fbeam)
  if (eerr) return
endif

do ihh = 1, run%n_harm
  first_was = ffield(ihh)%slip%first
  hsuf = ''
  if (ffield(ihh)%harm /= 1) write (hsuf, '(a, i0)') '-h', ffield(ihh)%harm

  if (run%ustate%active .or. .not. whole) then

    ! A scratch copy of the slices this frame writes, run%dump_is1 to run%dump_is2 of the
    ! window in time order. Each slice is taken from its record through fel_field_index.
    ! The copy is the slab and not the window. So a slice-range frame costs the slices it
    ! writes, and a slab across the record's rotation seam comes out in order. For a
    ! whole-window frame the copy is the window in time order, as the averaged path below
    ! rotates the live record to be.
    !
    ! Inside an unaveraged segment the record is the substep's midpoint field, half a
    ! substep past the plane this frame is named for (fel_unaveraged_mod). The frame
    ! carries the plane's own field. So the copy is carried back and written, and the live
    ! record is not touched. A run with frames on and one with them off then hold the same
    ! field to the bit. Carrying the live record back and forward again would leave two
    ! transform pairs of rounding in it, 5e-15 of the field on the mixed line
    ! (docs/validation.md). The carry is a transform pair a slice, with the phase table
    ! formed for each slice as fel_field_diffract_by forms it. The slices run in parallel as
    ! the step runs them, on the plans fel_field_kernel_init warmed at the segment's entry.
    ! Until 2026-10-30 the copy was the whole window, carried on one thread and then cut.
    ! That cost 2.5 s a frame for 12 slices of a 4212-slice window at 128 by 128
    ! (docs/performance.md).

    ! Under a field coarser than the slices (fel_field_sample) the slab is the planes the
    ! range's slices read, from the plane at or below the first to the one at or above
    ! the last, and its z origin is the first plane's.

    wf => ffield(ihh)%wf
    nslice = size(wf%Ex, 3)
    m = ffield(ihh)%n_per
    js1 = (run%dump_is1 - 1) / m + 1
    js2 = (run%dump_is2 - 1) / m + 1
    if (mod(run%dump_is2 - 1, m) /= 0) js2 = min(nslice, js2 + 1)
    nsub = js2 - js1 + 1
    if (allocated(wf_sub%Ex)) deallocate (wf_sub%Ex)
    if (allocated(wf_sub%Ey)) deallocate (wf_sub%Ey)
    allocate (wf_sub%Ex(size(wf%Ex, 1), size(wf%Ex, 2), nsub))
    if (allocated(wf%Ey)) allocate (wf_sub%Ey(size(wf%Ey, 1), size(wf%Ey, 2), nsub))
    wf_sub%dx = wf%dx
    wf_sub%dy = wf%dy
    wf_sub%dz = wf%dz
    wf_sub%wavelength = wf%wavelength
    wf_sub%s_position = wf%s_position

    !$OMP parallel do private(ifld)
    do j = 1, nsub
      ifld = fel_field_index(ffield(ihh)%slip, js1 + j - 1, nslice)
      wf_sub%Ex(:,:,j) = wf%Ex(:,:,ifld)
      if (allocated(wf%Ey)) wf_sub%Ey(:,:,j) = wf%Ey(:,:,ifld)
    enddo
    !$OMP end parallel do

    if (run%ustate%active) then
      any_err = .false.
      !$OMP parallel do private(eerr) reduction(.or.: any_err)
      do j = 1, nsub
        call fel_field_diffract_by (wf_sub, j, -0.5_rp * run%ustate%dsub, eerr)
        any_err = any_err .or. eerr
      enddo
      !$OMP end parallel do
      if (any_err) return
    endif

    ! The slab keeps its place. The mesh's z origin is the first slice's. So a reader
    ! that walks the standard attributes puts the cut field where the whole window's
    ! would sit. It sits beside the cut beam, whose timeOffset places each slice.

    call wavefront_write_openpmd (wf_sub, trim(prefix) // trim(hsuf) // '.wf.h5', eerr, &
                                  z_offset = (js1 - 1) * wf%dz)
    if (eerr) return
  else

    ! The whole window on the averaged path is written in place: rotated to time order,
    ! written, rotated back. A rotation is a permutation and exact. The field frame of a
    ! saturated run is the largest thing this program writes, and this path copies none
    ! of it.

    if (first_was /= 0) then
      ffield(ihh)%wf%Ex = cshift(ffield(ihh)%wf%Ex, shift = first_was, dim = 3)
      if (allocated(ffield(ihh)%wf%Ey)) &
          ffield(ihh)%wf%Ey = cshift(ffield(ihh)%wf%Ey, shift = first_was, dim = 3)
    endif
    call wavefront_write_openpmd (ffield(ihh)%wf, trim(prefix) // trim(hsuf) // '.wf.h5', eerr)
    if (first_was /= 0) then      ! Back to the rotation the walk carries, whatever the write did.
      ffield(ihh)%wf%Ex = cshift(ffield(ihh)%wf%Ex, shift = -first_was, dim = 3)
      if (allocated(ffield(ihh)%wf%Ey)) &
          ffield(ihh)%wf%Ey = cshift(ffield(ihh)%wf%Ey, shift = -first_was, dim = 3)
    endif
    if (eerr) return
  endif
  call fel_frame_attributes (trim(prefix) // trim(hsuf) // '.wf.h5', run, ie, eerr)
  if (eerr) return
  call fel_write_slippage_residual (trim(prefix) // trim(hsuf) // '.wf.h5', ffield(ihh)%slip%accuslip, eerr)
  if (eerr) return
  if (ckpt) then
    call fel_write_checkpoint_group (trim(prefix) // trim(hsuf) // '.wf.h5', run%lat%branch(0)%ele(ie)%s, ie, &
                                     run%lat%branch(0)%ele(ie)%name, eerr)
    if (eerr) return
  endif
enddo

err_flag = .false.

end subroutine fel_dump_frame

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_dump_field_set (run, prefix, err_flag, ix_done)
!
! Routine to write the whole field set at the given filename prefix as openPMD
! EXT_Wavefront (.wf.h5), both polarizations as components of one mesh record. Genesis
! .fld.h5 conversion lives in lucifer/tests/scripts/convert_genesis.py.
! The fundamental keeps the pre-harmonic names. A harmonic's files carry -h<h>.
! Every record is unrotated to time order first, since each field owns its rotation
! state. The fields move in lockstep, but the cshift must run per record.
!
! Input:
!   run       -- fel_run_struct: Run state with the field set to dump.
!   prefix    -- character(*): Filename prefix. Format-specific suffixes are appended.
!   ix_done   -- integer, optional: Index of the element the dump completes. The files
!                  carry the checkpoint group where that boundary is eligible
!                  (fel_checkpoint_eligible). If not present, no group is written.
!
! Output:
!   run       -- fel_run_struct: Field records unrotated to time order (slip%first = 0).
!   err_flag  -- logical: Set True if a file could not be written. False otherwise.
!-

subroutine fel_dump_field_set (run, prefix, err_flag, ix_done)

type (fel_run_struct), target :: run
type (fel_field_struct), pointer :: ffield(:)
type (ele_struct), pointer :: done
character(*) prefix
logical err_flag
integer ihh
integer, optional :: ix_done
logical eerr, ckpt
character(8) hsuf
real(rp), pointer :: z_now
integer n_harm
character(*), parameter :: r_name = 'fel_dump_field_set'

!

err_flag = .true.
ffield => run%ffield
z_now => run%z_now
n_harm = run%n_harm
ckpt = present(ix_done)
if (ckpt) ckpt = fel_checkpoint_eligible(run, ix_done)
if (ckpt) done => run%lat%branch(0)%ele(ix_done)

do ihh = 1, n_harm
  if (ffield(ihh)%slip%first /= 0) then
    ffield(ihh)%wf%Ex = cshift(ffield(ihh)%wf%Ex, shift = ffield(ihh)%slip%first, dim = 3)
    if (allocated(ffield(ihh)%wf%Ey)) &
        ffield(ihh)%wf%Ey = cshift(ffield(ihh)%wf%Ey, shift = ffield(ihh)%slip%first, dim = 3)
    ffield(ihh)%slip%first = 0
  endif

  hsuf = ''
  if (ffield(ihh)%harm /= 1) write (hsuf, '(a, i0)') '-h', ffield(ihh)%harm

  call wavefront_write_openpmd (ffield(ihh)%wf, prefix // trim(hsuf) // '.wf.h5', eerr)
  if (eerr) return
  call fel_write_slippage_residual (prefix // trim(hsuf) // '.wf.h5', ffield(ihh)%slip%accuslip, eerr)
  if (eerr) return
  if (ckpt) then
    call fel_write_checkpoint_group (prefix // trim(hsuf) // '.wf.h5', done%s, ix_done, done%name, eerr)
    if (eerr) return
  endif
enddo

err_flag = .false.

end subroutine fel_dump_field_set

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function si_text (value, unit) result (txt)
!
! Routine to write a value with its unit scaled to a prefix, for a box of the printout:
! 30.824 pC, 150.992 pm. The mantissa keeps at most three decimals and loses its trailing
! zeros down to one, the micro prefix is fel_micro$, which fel_print_line writes as the
! micro sign or as u, and zero is written as 0 and the unit.
!
! Input:
!   value -- real(rp): The value in the unit.
!   unit  -- character(*): The unit, m or C.
!
! Output:
!   txt   -- character(:), allocatable: The text.
!-

recursive function si_text (value, unit, sig) result (txt)

real(rp) value, mant
character(*) unit
integer, optional :: sig
character(:), allocatable :: txt, num, tail
character(14) str
character(12) short
integer ib

!

if (value == 0) then
  txt = '0 ' // unit
  return
endif

! A unit of two words, m rad, takes its prefix on the first.

ib = index(unit, ' ')
if (ib > 0) then
  txt = si_text(value, unit(1:ib-1), sig) // unit(ib:)
  return
endif

str = adjustl(fel_si_str(value, unit))
ib = index(trim(str), ' ')
if (ib == 0) then
  txt = trim(str)
  return
endif
num = str(1:ib-1)
tail = trim(adjustl(str(ib+1:)))
if (index(num, '.') > 0 .and. scan(num, 'Ee') == 0) then
  do
    if (num(len(num):len(num)) /= '0') exit
    if (num(len(num)-1:len(num)-1) == '.') exit
    num = num(1:len(num)-1)
  enddo
endif
if (len(tail) == len(unit) + 1) then
  if (tail(1:1) == 'u') tail = fel_micro$ // tail(2:)
endif

! With sig, the mantissa is rounded to four significant digits, for a column that is
! read at a glance.

if (present(sig) .and. scan(num, 'Ee') == 0) then
  read (num, *) mant
  if (abs(mant) < 10) then
    write (short, '(f0.3)') mant
  elseif (abs(mant) < 100) then
    write (short, '(f0.2)') mant
  else
    write (short, '(f0.1)') mant
  endif
  num = trim(short)
  if (num(1:1) == '.') num = '0' // num
endif
txt = num // ' ' // tail

end function si_text

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function number_text (x) result (txt)
!
! Routine to write a pure number compactly for a box of the printout: 0.19, 7.0, 22.3,
! 1253, 5.90e-4, 6.41e5.
!
! Input:
!   x   -- real(rp): The number.
!
! Output:
!   txt -- character(:), allocatable: The text.
!-

function number_text (x) result (txt)

real(rp) x, ax
character(:), allocatable :: txt
character(24) str
integer ie, iexp

!

ax = abs(x)
if (ax == 0) then
  txt = '0'
elseif (ax >= 1e4_rp .or. ax < 1e-2_rp) then
  write (str, '(es10.2e3)') x
  str = adjustl(str)
  ie = index(str, 'E')
  read (str(ie+1:), *) iexp
  txt = str(1:ie-1) // 'e' // int_str(iexp)
elseif (ax >= 100) then
  txt = int_str(nint(x))
elseif (ax >= 10) then
  write (str, '(f0.1)') x
  txt = trim(str)
else
  write (str, '(f0.2)') x
  txt = trim(str)
  if (txt(len(txt):len(txt)) == '0') txt = txt(1:len(txt)-1)
endif
if (txt(1:1) == '.') txt = '0' // txt
if (len(txt) > 1) then
  if (txt(1:2) == '-.') txt = '-0' // txt(2:)
endif

end function number_text

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function byte_text (bytes) result (txt)
!
! Routine to write a count of bytes in kB, MB or GB, a thousand to the step.
!
! Input:
!   bytes -- real(rp): The bytes.
!
! Output:
!   txt   -- character(:), allocatable: The text.
!-

function byte_text (bytes) result (txt)

real(rp) bytes
character(:), allocatable :: txt
character(24) str

!

if (bytes >= 1e9_rp) then
  write (str, '(f0.1, a)') bytes / 1e9_rp, ' GB'
elseif (bytes >= 1e6_rp) then
  write (str, '(f0.1, a)') bytes / 1e6_rp, ' MB'
else
  write (str, '(f0.1, a)') bytes / 1e3_rp, ' kB'
endif
txt = trim(str)
if (txt(1:1) == '.') txt = '0' // txt

end function byte_text

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_write_header (run)
!
! Routine to print the run's configuration as one framed block for a human to read:
! what lattice, what beam, what radiation, what is switched on, where the output goes.
! Display only. The same values, resolved, are echoed machine-readably into the stats
! file's meta/ group.
!
! Input:
!   run -- fel_run_struct: The fully set-up run state.
!-

subroutine fel_write_header (run)

type (fel_run_struct), target :: run
type (branch_struct), pointer :: branch
type (fel_key_struct), allocatable :: keys(:)
type (fel_column_struct) c_lat, c_win, c_beam, c_field, c_run, c_bmad, c_wake, c_sc, c_left, c_right
character(400) line
character(*), parameter :: r_name = 'lucifer'
character(8) date
character(10) time
character(64) host
character(:), allocatable :: stamp
real(rp) n_pstep, n_gpt, secs, q_max, cur_peak, l_1d, l_coop, window, slip, dx, half, eps, ds_first
real(rp) e_slice, e_beamlet, b_pulse, b_dev
integer nfel, n_omp, ie, i, nstep, n_particle, ngrid, nstep_unavg, n_grid(3), is, is_max, nbeamlet, n_rot
integer je_first, n_harm_lit
logical any_segment, any_storage, any_warning
character(16) backend

!

branch => run%lat%branch(0)
nfel = count(run%is_fel)
n_omp = 1
!$ n_omp = omp_get_max_threads()
call fel_schema_table (keys)

! The logo, with the versions, the time and the machine. The line with the versions
! moves with every commit, and the sample output on the user guide's page leaves it out.

call date_and_time (date, time)
stamp = date(1:4) // '-' // date(5:6) // '-' // date(7:8) // ' ' // time(1:2) // ':' // time(3:4) // ':' // time(5:6)
host = fel_print_host()
if (host /= '') stamp = stamp // ' ' // fel_dot$ // ' ' // trim(host)
call fel_print_logo ('FEL tracking in Bmad', &
      trim(lucifer_version$) // ' ' // fel_dot$ // ' Bmad ' // trim(bmad_version_date), stamp)
call out_io (s_blank$, r_name, '')

! The files, each across the box, since a path may be long.

call fel_box_rule (fel_rule_top$, 'files', '', .false., .false.)
call fel_box_full ('input     ' // trim(run%input_file))
call fel_box_full ('lattice   ' // trim(run%lat_file))
call fel_box_full ('resolved  ' // trim(run%global%out_root) // '.resolved.toml')

! The lattice's facts beside the window.

call fel_column_add (c_lat, 'length', si_text(branch%ele(branch%n_ele_track)%s, 'm'))
call fel_column_add (c_lat, 'elements', int_str(branch%n_ele_track))
call fel_column_add (c_lat, 'FEL segments', int_str(nfel))
write (line, '(f0.2)') run%gamma0
call fel_column_add (c_lat, 'gamma0', trim(line))

! The window as the run has it. A count the input file left to the program is marked.

call fel_column_add (c_win, 'n_wavelength', int_str(run%fbeam%n_wavelength), &
      origin_mark(key_index('window', 'n_wavelength'), .true.))
call fel_column_add (c_win, 'n_slice', int_str(run%nslice), origin_mark(key_index('window', 'n_slice'), run%nslice > 1))
call add_keys (c_win, 'window', '', skip = [character(32) :: 'n_wavelength', 'n_slice'])
if (run%nslice > 1 .and. run%ffield(1)%n_per > 1) then
  call fel_column_add (c_win, 'field planes', int_str(run%nfield))
endif

call fel_box_rule (fel_rule_mid$, 'lattice', 'window', .false., .true.)
call fel_box_columns (c_lat, c_win)

! An FEL segment with an offset, a roll-off wavenumber or a kick says so, since each is
! read from the lattice and changes the tracking (fel-physics.md sec-element).

any_segment = .false.
do ie = 1, branch%n_ele_track
  if (.not. run%is_fel(ie)) cycle
  if (.not. run%und_of(ie)%field_general) cycle
  if (.not. any_segment) call fel_box_rule (fel_rule_mid$, 'segments', '', .true., .false.)
  any_segment = .true.
  write (line, '(a16, 4a)') branch%ele(ie)%name, 'x_offset ', si_text(run%und_of(ie)%ax, 'm'), &
        '   y_offset ', si_text(run%und_of(ie)%ay, 'm')
  call fel_box_full (trim(line))
  write (line, '(16x, a, f0.3, 5a)') 'kx ', run%und_of(ie)%kxb, ' /m   hkick ', &
        si_text(run%und_of(ie)%hkick, 'rad'), '   vkick ', si_text(run%und_of(ie)%vkick, 'rad')
  call fel_box_full (trim(line))
enddo

! The beam beside the field, with the keys of the noise in a part of their own under
! the field's.

call add_keys (c_beam, 'beam', '', skip = noise_keys)
call add_keys (c_beam, 'beam.resample', 'resample.')
call add_keys (c_field, 'field', '')
call fel_column_rule (c_field, 'noise')
call add_keys (c_field, 'beam', '', only = noise_keys)

call fel_box_rule (fel_rule_mid$, 'beam', 'field', .not. any_segment, .true.)
call fel_box_columns (c_beam, c_field)

call add_keys (c_run, 'run', '')
call add_keys (c_bmad, 'bmad', '')
call fel_box_rule (fel_rule_mid$, 'run', 'bmad', .true., .true.)
call fel_box_columns (c_run, c_bmad)

call add_keys (c_wake, 'wake', '')
call add_keys (c_sc, 'space_charge', '')
if (c_wake%n > 1 .or. c_sc%n > 1 .or. run%coll%wake%on) then
  call fel_box_rule (fel_rule_mid$, 'wake', 'space_charge', .true., .true.)
  call fel_box_columns (c_wake, c_sc)
endif

! What follows from the settings: the quantities a reader would otherwise work out by
! hand, with the ratios that say whether a setting is in proportion on the right.

n_grid = wavefront_shape(run%ffield(1)%wf)
dx = run%ffield(1)%wf%dx
half = 0.5_rp * dx * (n_grid(1) - 1)

q_max = 0;  is_max = 1
do is = 1, run%nslice
  if (sum(run%fbeam%slice(is)%weight(1:run%fbeam%slice(is)%n)) > q_max) then
    q_max = sum(run%fbeam%slice(is)%weight(1:run%fbeam%slice(is)%n))
    is_max = is
  endif
enddo
cur_peak = 0
if (run%fbeam%slice_spacing > 0) cur_peak = q_max * c_light / run%fbeam%slice_spacing
e_slice = q_max / e_charge
nbeamlet = run%fbeam%slice(is_max)%n / max(1, run%bparam%beamlet_size)

je_first = 0
do ie = 1, branch%n_ele_track
  if (run%is_fel(ie)) then
    je_first = ie
    exit
  endif
enddo
l_1d = 0;  l_coop = 0;  ds_first = 0
if (run%rho_1d > 0 .and. je_first > 0) then
  l_1d = twopi / (run%und_of(je_first)%ku * 4 * pi * sqrt(3.0_rp) * run%rho_1d)
  l_coop = run%winit%lambda0 / (4 * pi * run%rho_1d)
  ds_first = branch%ele(je_first)%value(ds_step$)
endif
window = run%nslice * run%fbeam%slice_spacing
slip = 0
if (run%nslice > 1 .and. allocated(run%ele_slip)) slip = (run%slip_und + sum(run%ele_slip)) * run%winit%lambda0
eps = 0
if (run%gamma0 > 0) eps = 0.5_rp * (run%beam_init%a_norm_emit + run%beam_init%b_norm_emit) / run%gamma0

if (run%rho_1d > 0) then
  call fel_column_add (c_left, 'rho (1D)', number_text(run%rho_1d))
  call fel_column_add (c_right, 'gain length (1D)', si_text(l_1d, 'm', 4))
endif
call fel_column_add (c_left, 'peak current', si_text(cur_peak, 'A', 4))
if (l_coop > 0) call fel_column_add (c_right, 'cooperation length', si_text(l_coop, 'm', 4))
if (run%sigma_beam > 0) then
  call fel_column_add (c_left, 'beam size', si_text(run%sigma_beam, 'm', 4))
  if (eps > 0 .and. run%winit%lambda0 > 0) call fel_column_add (c_right, '4 pi emittance / lambda0', &
        number_text(4 * pi * eps / run%winit%lambda0))
endif
if (run%nslice > 1) then
  call fel_column_add (c_left, 'window', si_text(window, 'm', 4))
  call fel_column_add (c_right, 'slippage / window', number_text(slip / window))
  call fel_column_add (c_left, 'slice spacing', si_text(run%fbeam%slice_spacing, 'm', 4))
  if (l_coop > 0) call fel_column_add (c_right, 'cooperation length / slice', &
        number_text(l_coop / run%fbeam%slice_spacing))
endif
call fel_column_add (c_left, 'grid cell', si_text(dx, 'm', 4))
if (run%sigma_beam > 0) then
  call fel_column_add (c_right, 'beam size / grid cell', number_text(run%sigma_beam / dx))
  call fel_column_add (c_right, 'grid half width / beam size', number_text(half / run%sigma_beam))
endif
if (run%global%source_filter .and. run%filter_edge > 0) then
  call fel_column_add (c_left, 'filter edge', si_text(run%filter_edge, 'rad', 4))
endif
call fel_column_add (c_left, 'electrons / slice', number_text(e_slice))
call fel_cost_count (run, nstep, n_particle, n_pstep, n_gpt, ngrid, nstep_unavg)
call fel_column_add (c_left, 'macroparticles', int_str(n_particle))
if (run%bparam%quiet_start .and. run%bparam%beam_file == '' .and. nbeamlet > 0) then
  call fel_column_add (c_left, 'beamlets / slice', int_str(nbeamlet))
  e_beamlet = e_slice / nbeamlet
  call fel_column_add (c_right, 'electrons / beamlet', number_text(e_beamlet))
endif
if (l_1d > 0 .and. ds_first > 0) call fel_column_add (c_right, 'gain length / step', number_text(l_1d / ds_first))

call fel_box_rule (fel_rule_mid$, 'derived', '', .true., .true.)
call fel_box_columns (c_left, c_right)

! The same numbers as values, for the stats file (fel_record_mod).

call fel_record_real ('setup', 'rho_1d', '1', 'The Pierce parameter of the 1D theory, at the peak current.', run%rho_1d)
call fel_record_real ('setup', 'gain_length_1d', 'm', 'The power gain length of the 1D theory, ' // &
      'lambda_u / (4 pi sqrt(3) rho), in the first FEL element.', l_1d)
call fel_record_real ('setup', 'cooperation_length', 'm', 'The cooperation length, lambda0 / (4 pi rho).', l_coop)
call fel_record_real ('setup', 'peak_current', 'A', 'The current of the slice with the most charge.', cur_peak)
call fel_record_real ('setup', 'beam_size', 'm', 'The rms beam size the grid was derived from.', run%sigma_beam)
call fel_record_real ('setup', 'window', 'm', 'The length of the time window, n_slice times the slice spacing.', window)
call fel_record_real ('setup', 'slippage', 'm', 'The slippage of the light over the line.', slip)
call fel_record_real ('setup', 'slippage_undulator', '1', &
      'The slippage in undulator steps (wavelengths).', run%slip_und)
if (allocated(run%ele_slip)) call fel_record_real ('setup', 'slippage_wavelengths', '1', &
      'The slippage of the light over the line (wavelengths), the sum the schedule applies.', &
      run%slip_und + sum(run%ele_slip))
call fel_record_int ('setup', 'slippage_interlude', 'The slippage in interludes and at the end (wavelengths).', &
      run%slip_int)
call fel_record_int ('setup', 'slippage_chicane', 'The slippage in chicanes (wavelengths).', run%slip_chic)
call fel_record_real ('setup', 'headroom_wavelengths', '1', 'The slippage a derived window holds headroom ' // &
      'for (wavelengths), 0 for a window the input set.', run%slip_headroom)
call fel_record_int ('setup', 'n_slice_headroom', 'The slices of headroom in a derived window.', run%n_headroom)
call fel_record_int ('setup', 'n_field_plane', 'The planes of the field over the window.', run%nfield)
call fel_record_int ('setup', 'field_n_wavelength', 'The wavelengths one plane of the field spans.', &
      run%ffield(1)%n_wavelength)
call fel_record_int ('setup', 'n_pad_free', 'The free slots of the spectral propagator''s pad along tau.', &
      run%spectral%n_free)
call fel_record_int ('setup', 'n_pad_sponge', 'The sponge slots of the spectral propagator''s pad along tau.', &
      run%spectral%n_sponge)
call fel_record_int ('setup', 'slices_per_plane', 'The slices one plane of the field spans.', run%ffield(1)%n_per)
call fel_record_int ('setup', 'n_interlude_piece', 'The pieces run.interlude_ds_step cuts interludes into, ' // &
      '0 where every interlude is tracked whole.', merge(run%n_int_piece, 0, run%n_int_ele > 0))
call fel_record_int ('setup', 'n_interlude_element', 'The elements run.interlude_ds_step cuts.', run%n_int_ele)
call fel_record_real ('setup', 'grid_cell', 'm', 'The transverse grid''s cell.', dx)
call fel_record_real ('setup', 'grid_half_width', 'm', 'The transverse grid''s half width.', half)
call fel_record_int ('setup', 'grid_n_pts', 'The transverse grid''s points on a side.', n_grid(1))
call fel_record_real ('setup', 'filter_edge', 'rad', 'The source filter''s edge angle, 0 with the filter off.', &
      merge(run%filter_edge, 0.0_rp, run%global%source_filter))
call fel_record_real ('setup', 'electrons_per_slice', '1', 'The electrons in the slice with the most charge.', e_slice)
call fel_record_int ('setup', 'n_macroparticle', 'The macroparticles of the whole bunch.', n_particle)
call fel_record_int ('setup', 'beamlets_per_slice', 'The beamlets in the slice with the most charge.', nbeamlet)
call fel_record_int ('setup', 'beamlet_size', 'The macroparticles of one beamlet.', run%bparam%beamlet_size)
call fel_record_real ('setup', 'step', 'm', 'The integration step of the first FEL element.', ds_first)

! What the run will write and hold, stated before it tracks: the pulse, the frames, the
! field's reductions and the device's resident bytes.

call fel_cost_count (run, nstep, n_particle, n_pstep, n_gpt, ngrid, nstep_unavg)
call bank_bytes (n_rot, b_pulse, n_harm_lit)
any_storage = .false.
if (n_harm_lit > 0) then
  call storage_rule ()
  if (run%global%write_escaped) then
    write (line, '(3a, i0, a)') 'pulse and escaped files   ', trim(byte_text(b_pulse)), '   ', n_rot, &
          ' slices leave the window'
  else
    write (line, '(3a, i0, a)') 'pulse file   ', trim(byte_text(b_pulse)), '   ', n_rot, ' slices leave the window'
  endif
  call fel_box_full (trim(line))
endif
if (run%global%dump_at_comb) then
  call storage_rule ()
  write (line, '(3a, i0, a, i0, a)') 'frames   ', trim(byte_text(run%stats%nrec * frame_bytes())), '   ', &
        run%stats%nrec, ' at the comb over ', run%dump_is2 - run%dump_is1 + 1, ' slices'
  call fel_box_full (trim(line))
endif
if (run%stats%n_red > 0) then
  call storage_rule ()
  write (line, '(3a, i0, a)') 'reductions   ', trim(byte_text(8.0_rp * run%stats%nrec * run%stats%n_red * &
        (real(ngrid, rp)**2 + 2.0_rp * run%nslice * ngrid))), '   ', run%stats%n_red, ' of the field set, in the stats file'
  call fel_box_full (trim(line))
endif
if (run%dev%on) then
  call storage_rule ()
  call device_bytes (b_dev)
  write (line, '(5a, i0, a)') 'device   ', trim(byte_text(b_dev)), ' resident on ', trim(run%dev%name), '   ', &
        run%dev%nplane, ' planes'
  call fel_box_full (trim(line))
  write (line, '(a, es8.2, a, es8.2, a)') 'device   FP32 transform pair loss ', fel_dev_pair_loss(ngrid), &
        ' a pair, ', fel_dev_pair_loss(ngrid) * nstep, ' over the line'
  call fel_box_full (trim(line))
endif

! The work: the counts are exact, and the estimate is one machine's, which the row names
! (fel_cost_mod). load_only stops after this block, so an input file's work can be seen
! without tracking it.

call fel_cost_estimate (run, secs, backend)
call fel_record_int ('setup', 'n_slice_leaving', 'The slices of light that leave the window over the line.', n_rot)
call fel_record_real ('setup', 'pulse_bytes', 'B', 'The size of the pulse files, as stated before tracking.', &
      merge(b_pulse, 0.0_rp, n_harm_lit > 0))
if (run%global%dump_at_comb) then
  call fel_record_int ('setup', 'n_frame', 'The frames the run writes at the comb.', run%stats%nrec)
  call fel_record_real ('setup', 'frame_bytes', 'B', 'The size of all frames together.', run%stats%nrec * frame_bytes())
endif
call fel_record_int ('timing', 'n_thread', 'The OpenMP threads of the run.', n_omp)
call fel_record_int ('setup', 'n_fel_step', 'The FEL steps of the line.', nstep)
call fel_record_int ('setup', 'n_unaveraged_substep', 'The substeps of the unaveraged elements.', nstep_unavg)
if (run%dev%on) call fel_record_real ('setup', 'device_pair_loss', '1', 'The projected loss of field energy ' // &
      'in one FP32 transform pair on the device.', fel_dev_pair_loss(ngrid))
call fel_record_real ('setup', 'particle_steps', '1', 'Macroparticles times FEL steps.', n_pstep)
call fel_record_real ('timing', 'estimate_walk', 's', 'The estimate of the walk''s time, on the machine ' // &
      'estimate_machine names.', secs)
call fel_record_str ('timing', 'estimate_path', 'The path the estimate is for, CPU or device.', trim(backend))
call fel_record_str ('timing', 'estimate_machine', 'The machine the estimate''s rates were timed on.', fel_cost_machine$)
call fel_box_rule (fel_rule_mid$, 'work', '', .not. any_storage, .false.)
write (line, '(a, i0, a, i0, 2a)') 'threads ', n_omp, '      FEL steps ', nstep, &
      '      particle-steps ', number_text(n_pstep)
call fel_box_full (trim(line))
call fel_box_full ('estimate ' // trim(fel_cost_secs_str(secs)) // ', ' // trim(backend) // &
      ' path, as timed on ' // fel_cost_machine$)
call fel_box_rule (fel_rule_bottom$, '', '', .false., .false.)
call fel_print_line ('   * derived     ~ command line     (blank) input file or default')

call write_warnings ()
call out_io (s_blank$, r_name, '')

!-------------------------------------------------------
contains

!+
! Function key_index (table, name) result (ik)
!
! Routine to return a key's place in the schema's table, 0 where there is no such key.
!-

function key_index (table, name) result (ik)

character(*) table, name
integer ik

!

do ik = 1, size(keys)
  if (keys(ik)%table == table .and. keys(ik)%name == name) return
enddo
ik = 0

end function key_index

!------------------------------------------------------------------------------
!+
! Function origin_mark (ik, computed) result (mark)
!
! Routine to return the provenance mark of key ik: * where the program derived its value,
! ~ where the command line set it, and blank for the input file or a default. A key the
! input left unset whose value the run computed, the slices of a derived window among
! them, is marked derived when computed is true.
!-

function origin_mark (ik, computed) result (mark)

integer ik
logical computed
character(1) mark

!

mark = ' '
if (ik == 0 .or. .not. allocated(run%key_origin)) return
select case (run%key_origin(ik))
case (fel_origin_derived$);  mark = '*'
case (fel_origin_command$);  mark = '~'
case (fel_origin_default$);  if (computed) mark = '*'
end select

end function origin_mark

!------------------------------------------------------------------------------
!+
! Subroutine add_keys (col, table, prefix, skip, only)
!
! Routine to add the keys of one table to a column: every key the input file or the
! command line set or the program derived, and the defaults of shown_defaults, which
! decide how a run behaves. The value is written by key_text, from the key's own kind
! and unit.
!
! Input:
!   col    -- fel_column_struct: The column.
!   table  -- character(*): The table.
!   prefix -- character(*): Written before each key's name, resample. for the keys of
!               [beam.resample].
!   skip   -- character(32)(:), optional: Keys to leave out.
!   only   -- character(32)(:), optional: The only keys to take.
!-

subroutine add_keys (col, table, prefix, skip, only)

type (fel_column_struct) col
character(*) table, prefix
character(32), optional :: skip(:), only(:)
integer ik
character(48) dotted

!

do ik = 1, size(keys)
  if (keys(ik)%table /= table) cycle
  if (present(skip)) then
    if (any(skip == keys(ik)%name)) cycle
  endif
  if (present(only)) then
    if (.not. any(only == keys(ik)%name)) cycle
  endif
  dotted = trim(table) // '.' // trim(keys(ik)%name)
  if (run%key_origin(ik) == fel_origin_default$ .and. .not. any(shown_defaults == dotted)) cycle
  call fel_column_add (col, prefix // trim(keys(ik)%name), key_text(keys(ik)), origin_mark(ik, .false.))
enddo

end subroutine add_keys

!------------------------------------------------------------------------------
!+
! Function key_text (key) result (txt)
!
! Routine to write a key's value for the header from its kind and its unit: a real with
! a unit that takes a prefix scaled to it, a string without its quotes, a list with
! commas, and anything else as the input file writes it.
!-

function key_text (key) result (txt)

type (fel_key_struct) key
character(:), allocatable :: txt, sv
real(rp) rv
integer iv, j
logical lv
integer, allocatable :: ilist(:)
character(400), allocatable :: slist(:)

!

call fel_schema_value (run, key, rv, iv, lv, sv, ilist, slist)

select case (key%kind)
case (key_real$)
  if (any(scaled_units == key%unit)) then
    txt = si_text(rv, trim(key%unit))
  else
    txt = fel_schema_text(run, key)
    if (key%unit /= '') txt = txt // ' ' // trim(key%unit)
  endif
case (key_string$, key_path$)
  txt = sv
case (key_int_list$)
  txt = ''
  do j = 1, size(ilist)
    if (j > 1) txt = txt // ', '
    txt = txt // int_str(ilist(j))
  enddo
case (key_string_list$)
  txt = ''
  do j = 1, size(slist)
    if (j > 1) txt = txt // ', '
    txt = txt // trim(slist(j))
  enddo
case default
  txt = fel_schema_text(run, key)
end select

end function key_text

!------------------------------------------------------------------------------
!+
! Subroutine storage_rule ()
!
! Routine to write the rule above the first row of what the run stores, once.
!-

subroutine storage_rule ()

if (any_storage) return
call fel_box_rule (fel_rule_mid$, 'storage', '', .true., .false.)
any_storage = .true.

end subroutine storage_rule

!------------------------------------------------------------------------------
!+
! Function frame_bytes () result (bytes)
!
! Routine to return the bytes of one frame of the series: 68 a macroparticle and 16 a
! grid point, over the slice range the frames keep.
!-

function frame_bytes () result (bytes)

real(rp) bytes, n_rng

!

n_rng = real(run%dump_is2 - run%dump_is1 + 1, rp) / max(1, run%nslice)
bytes = (68.0_rp * n_particle + 16.0_rp * n_gpt) * n_rng

end function frame_bytes

!------------------------------------------------------------------------------
!+
! Subroutine bank_bytes (n_rot, bytes, n_lit)
!
! Routine to return what the escaped-field bank will write: the slices that leave the
! window over the tracked line, the bytes of the files the run keeps over the members of
! the field set that carry light, and how many members those are. The banked planes
! stream into those files as they leave (fel_drain_bank), so the bytes kept are also the
! peak on disk.
!
! The slices are the rotations the walk makes over the tracked elements, from the
! schedule's own slippage (fel_setup_slippage) and the residual a continuation starts
! with: the record rotates once |accuslip| passes 0.8 of a slice, so S wavelengths give
! the least integer at or above S/n_wavelength - 0.8. Light leaves only if one of the
! last n_rot slices holds charge or field at the start. A file holds 16 bytes a grid
! point a component, a plane per slice, and a row of 25 params a component and the
! transmission position per banked slice. The pulse holds the window's slices too.
!-

subroutine bank_bytes (n_rot, bytes, n_lit)

real(rp) bytes, s_line, b_plane, b_rows
integer n_rot, n_lit, ihh, is1, if1, js, ng(3), n_comp
logical lit

!

n_rot = 0;  bytes = 0;  n_lit = 0
if (.not. run%global%keep_escaped_field .or. run%nslice <= 1) return
if (.not. allocated(run%und_slip)) return
s_line = run%ffield(1)%slip%accuslip + sum(run%und_slip(run%i_start:run%i_end)) + &
         sum(run%ele_slip(run%i_start:run%i_end))
n_rot = max(0, ceiling(s_line / run%fbeam%n_wavelength - 0.8_rp))

! The spectral propagator banks a plane at each whole plane the ramp crosses
! (fel_spectral_apply), the field's own spacing, and its count starts from the ramp's
! own residual.

if (run%spectral%on) then
  s_line = run%ffield(1)%slip%ramp_slip + sum(run%und_slip(run%i_start:run%i_end)) + &
           sum(run%ele_slip(run%i_start:run%i_end))
  n_rot = max(0, floor(s_line / run%ffield(1)%n_wavelength + 1e-9_rp))
endif
if (n_rot == 0) return
is1 = max(1, run%nslice - n_rot * run%ffield(1)%n_per + 1)
if1 = max(1, run%nfield - n_rot + 1)

do ihh = 1, run%n_harm
  lit = .false.
  do js = is1, run%nslice
    if (sum(run%fbeam%slice(js)%weight(1:run%fbeam%slice(js)%n)) > 0) lit = .true.
  enddo
  if (any(run%ffield(ihh)%wf%Ex(:,:,if1:run%nfield) /= 0)) lit = .true.
  if (allocated(run%ffield(ihh)%wf%Ey)) then
    if (any(run%ffield(ihh)%wf%Ey(:,:,if1:run%nfield) /= 0)) lit = .true.
  endif
  if (.not. lit) cycle
  n_lit = n_lit + 1

  ng = wavefront_shape(run%ffield(ihh)%wf)
  n_comp = 1
  if (allocated(run%ffield(ihh)%wf%Ey)) n_comp = 2
  b_plane = 16.0_rp * ng(1) * ng(2) * n_comp
  b_rows = 8.0_rp * n_rot * (1 + 25 * n_comp)
  bytes = bytes + (run%nfield + n_rot) * b_plane + b_rows
  if (run%global%write_escaped) bytes = bytes + n_rot * b_plane + b_rows
  call fel_record_real ('setup', 'bank_bytes_h' // int_str(run%ffield(ihh)%harm), 'B', &
        'The size of this harmonic''s pulse file, and of its escaped file with it under ' // &
        'run.write_escaped, as stated before tracking.', (run%nfield + n_rot) * b_plane + b_rows + &
        merge(n_rot * b_plane + b_rows, 0.0_rp, run%global%write_escaped))
enddo

end subroutine bank_bytes

!------------------------------------------------------------------------------
!+
! Subroutine device_bytes (bytes)
!
! Routine to return the device's resident bytes (fel_device_resident_bytes).
!-

subroutine device_bytes (bytes)

use fel_device_mod, only: fel_device_resident_bytes

real(rp) bytes

!

bytes = fel_device_resident_bytes(run%dev)

end subroutine device_bytes

!------------------------------------------------------------------------------
!+
! Subroutine write_warnings ()
!
! Routine to write, under the header, a line for each quantity of the derived part that
! is out of proportion, each beginning with an exclamation mark so that a log can be
! searched for them. Each bound is half of what the setup derives for itself where the
! input file states nothing, or a count the loader already acts on.
!-

subroutine write_warnings ()

!

any_warning = .false.
if (run%nslice > 1 .and. window > 0) then
  if (slip / window > 0.5_rp) call warn ('slippage / window = ' // number_text(slip / window) // &
        ', above 0.5')
  if (l_coop > 0) then
    if (l_coop / run%fbeam%slice_spacing < 2) call warn ('cooperation length / slice = ' // &
          number_text(l_coop / run%fbeam%slice_spacing) // ', under 2')
  endif
endif
if (run%sigma_beam > 0) then
  if (run%sigma_beam / dx < 3.5_rp) call warn ('beam size / grid cell = ' // &
        number_text(run%sigma_beam / dx) // ', under 3.5')
  if (half / run%sigma_beam < 4.5_rp) call warn ('grid half width / beam size = ' // &
        number_text(half / run%sigma_beam) // ', under 4.5')
endif
if (l_1d > 0 .and. ds_first > 0) then
  if (l_1d / ds_first < 5) call warn ('gain length / step = ' // number_text(l_1d / ds_first) // ', under 5')
endif
if (run%bparam%quiet_start .and. run%bparam%beam_file == '' .and. nbeamlet > 0) then
  if (nbeamlet < 128) call warn ('beamlets / slice = ' // int_str(nbeamlet) // ', under 128')
  if (run%bparam%shot_noise .and. e_slice / nbeamlet < 1) call warn ('electrons / beamlet = ' // &
        number_text(e_slice / nbeamlet) // ', under 1')
endif
if (run%n_harm > 1 .and. .not. run%bparam%shot_noise .and. run%bparam%quiet_start) then
  call warn ('a harmonic is carried with shot_noise = false')
endif

end subroutine write_warnings

!------------------------------------------------------------------------------
!+
! Subroutine warn (text)
!
! Routine to write one warning line, with a blank line before the first.
!-

subroutine warn (text)

character(*) text

!

if (.not. any_warning) call out_io (s_blank$, r_name, '')
any_warning = .true.
call out_io (s_blank$, r_name, ' ! ' // text)

end subroutine warn

end subroutine fel_write_header

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_write_footer (run)
!
! Routine to close a run for a human: how far it went and how long it took, what the
! light ended up as, and which files were written. The file list comes last, so the last
! thing on screen is the list of what to open.
!
! The list is also what lucifer.runner.Result reads to find the run's dumps, so it names
! every dump this run wrote and no other: the initial state, the dumps at elements, the
! final state, the pulse and, with write_escaped, the escaped file, and the frames. A
! frame series is one
! line per kind of file, its record numbers written as a shell brace range,
! <out_root>-{000001..000012}.beam.h5, with the size of the whole series.
!
! Input:
!   run -- fel_run_struct: The finished run state.
!-

subroutine fel_write_footer (run)

type (fel_run_struct), target :: run
type (fel_stats_struct), pointer :: stats
type (fel_convergence_struct) cvg
type (fel_column_struct) c_left, c_right
character(400) line
character(*), parameter :: r_name = 'lucifer'
integer ir, ih, n_listed, is, n_thin, n_shown, nbl, ie, nrec, n_left, ip
real(rp) pow, ene, bun, pul, worst, secs, ratio_s, e_thin, walk, t_phase(8), t_other
logical few, any_warning, rec_err
character(8) hsuf
character(16) backend
character(420) root, dname
character(20), parameter :: phase_name(8) = [character(20) :: 'particle push', 'field solve', &
      'unaveraged step', 'device step', 'slippage', 'interludes', 'stats and diag', 'dumps']
integer, parameter :: n_thin_max$ = 8

!

stats => run%stats
call out_io (s_blank$, r_name, '')

! The light at the exit, as fel_stats_exit_light takes it, and the split of the power
! at the mode's angle where the mode peaked and at the last record
! (fel_stats_convergence).

call fel_stats_exit_light (stats, pow, ene, bun, pul)
call fel_stats_convergence (stats, cvg)

call fel_column_add (c_left, 'exit power', si_text(pow, 'W'))
call fel_column_add (c_right, 'exit pulse energy', si_text(pul, 'J'))
write (line, '(f6.4)') bun
call fel_column_add (c_left, 'exit bunching', trim(adjustl(line)))
call fel_column_add (c_right, 'window energy', si_text(ene, 'J'))
write (line, '(f0.3, a)') run%z_now, ' m'
call fel_column_add (c_left, 'tracked to', trim(line))
call fel_column_add (c_right, 'element ends', int_str(run%stats%iend))
if (cvg%ok) then
  call fel_column_add (c_left, 'split angle', si_text(cvg%angle, 'rad'))
  write (line, '(f0.3, a)') cvg%z(1), ' m'
  call fel_column_add (c_right, 'mode peak at', trim(line))
  call fel_column_add (c_left, 'outside / inside', number_text(cvg%ratio(1)) // ' (peak)')
  call fel_column_add (c_right, 'outside / inside', number_text(cvg%ratio(2)) // ' (exit)')
endif
if (run%ffield(1)%slip%timerun) then
  do ih = 1, run%n_harm
    write (line, '(a, i0)') 'escaped energy h', run%ffield(ih)%harm
    call fel_column_add (c_left, trim(line), si_text(run%ffield(ih)%slip%u_escaped, 'J'))
    call fel_column_add (c_right, 'slices that left', int_str(run%n_banked(ih)))
  enddo
endif
if (run%global%migrate .and. run%nslice > 1) then
  call fel_column_add (c_left, 'particles migrated', int_str(run%n_moved_tot))
  call fel_column_add (c_right, 'charge dropped', si_text(run%charge_dropped_tot, 'C'))
endif
if (run%fp32%on) then
  call fel_column_add (c_left, 'FP32 worst pz', number_text(run%fp32%worst(5)))
  call fel_column_add (c_right, 'FP32 worst theta', number_text(run%fp32%worst(6)))
  call fel_column_add (c_left, 'FP32 worst phasor', number_text(run%fp32%worst(7)))
  call fel_column_add (c_right, 'FP32 guard ulp', number_text(run%fp32%ulp_min))
endif

call fel_box_rule (fel_rule_top$, 'result', '', .false., .true.)
call fel_box_columns (c_left, c_right)

call fel_record_real ('run', 'exit_power', 'W', 'The peak power in the window at the last record.', pow)
call fel_record_real ('run', 'exit_window_energy', 'J', 'The energy in the window at the last record.', ene)
call fel_record_real ('run', 'exit_pulse_energy', 'J', 'The energy of the pulse, the window with the light ' // &
      'that left it.', pul)
call fel_record_real ('run', 'exit_bunching', '1', 'The bunching at the last record.', bun)
call fel_record_real ('run', 's_end', 'm', 'Where tracking ended.', run%z_now)
if (run%ffield(1)%slip%timerun) then
  do ih = 1, run%n_harm
    call fel_record_real ('run', 'escaped_energy_h' // int_str(run%ffield(ih)%harm), 'J', &
          'The energy of this harmonic that slippage carried out of the window over the run.', &
          run%ffield(ih)%slip%u_escaped)
  enddo
endif

! The time, by the phases that took any, with the walk against its estimate. --verbose
! prints the table of every phase as well.

walk = fel_timer_seconds(fel_t_walk$)
t_phase(1) = fel_timer_seconds(fel_t_particles$)
t_phase(2) = fel_timer_seconds(fel_t_solve$) + fel_timer_seconds(fel_t_source$)
t_phase(3) = fel_timer_seconds(fel_t_unavg$)
t_phase(4) = fel_timer_seconds(fel_t_device$)
t_phase(5) = fel_timer_seconds(fel_t_slippage$)
t_phase(6) = fel_timer_seconds(fel_t_seam$) + fel_timer_seconds(fel_t_drift$) + fel_timer_seconds(fel_t_gen_int$)
t_phase(7) = fel_timer_seconds(fel_t_stats$) + fel_timer_seconds(fel_t_ele_end$)
t_phase(8) = fel_timer_seconds(fel_t_dumps$)
t_other = walk - sum(t_phase)

c_left%n = 0;  c_right%n = 0
n_left = 0
do ip = 1, size(t_phase)
  if (walk <= 0 .or. t_phase(ip) < 0.005_rp * walk) then
    t_other = t_other + t_phase(ip)
    cycle
  endif
  call add_time (phase_name(ip), t_phase(ip))
enddo
! Every phase of the timer with the times it was entered, for the stats file. The boxes
! above merge some of them.

do ip = 1, fel_t_n$
  call fel_record_real ('timing', trim(record_name(fel_phase_name(ip))), 's', &
        'Wall-clock time of the walk in this phase.', fel_timer_seconds(ip))
  call fel_record_int ('timing', trim(record_name(fel_phase_name(ip))) // '_regions', &
        'The times this phase was entered.', fel_timer_regions(ip))
enddo
call fel_record_real ('timing', 'walk', 's', 'Wall-clock time of the walk.', walk)
call fel_record_real ('timing', 'run', 's', 'Wall-clock time of the run up to the closing printout.', &
      fel_timer_run_seconds())
if (walk > 0 .and. t_other >= 0.005_rp * walk) call add_time ('other', t_other)

call fel_cost_estimate (run, secs, backend)
if (mod(n_left, 2) == 1) call fel_column_add (c_right, '', '')
call fel_column_add (c_left, 'walk', secs_right(walk))
if (secs > 0) then
  call fel_column_add (c_right, 'walk / estimate', number_text(walk / secs))
else
  call fel_column_add (c_right, '', '')
endif
call fel_column_add (c_left, 'run', secs_right(fel_timer_run_seconds()))

call fel_box_rule (fel_rule_mid$, 'timing', '', .true., .true.)
call fel_box_columns (c_left, c_right)
if (fel_print_com%verbose) call fel_timer_write (r_name)

! Every file the run wrote. A dump is gated by the switch that writes it and then by its
! existence, so a file from an earlier run with the same out_root is not listed. A frame
! series is one row in shell brace syntax.

call fel_box_rule (fel_rule_mid$, 'output', '', .true., .false.)
n_listed = 0
root = run%global%out_root
if (run%global%write_initial) call note_dumps (trim(root) // '-initial', .true., .true.)
do ie = run%i_start, run%i_end
  if (.not. (run%dump_beam_here(ie) .or. run%dump_field_here(ie))) cycle
  write (dname, '(2a, i0, 2a)') trim(root), '-at', ie, '-', trim(run%lat%branch(0)%ele(ie)%name)
  call note_dumps (trim(dname), run%dump_beam_here(ie), run%dump_field_here(ie))
enddo
call note_dumps (trim(root) // '-final', .true., .true.)
do ih = 1, run%n_harm
  if (run%n_banked(ih) == 0) cycle
  if (run%global%write_escaped) call note_file (trim(fel_bank_file_name(run, ih, 'escaped')))
  call note_file (trim(fel_bank_file_name(run, ih, 'pulse')))
enddo
nrec = run%stats%irec
if (run%global%dump_at_comb .and. nrec > 0) then
  call note_series ('.beam.h5')
  do ih = 1, run%n_harm
    hsuf = ''
    if (run%ffield(ih)%harm /= 1) write (hsuf, '(a, i0)') '-h', run%ffield(ih)%harm
    call note_series (trim(hsuf) // '.wf.h5')
  enddo
endif
call note_file (trim(run%global%out_root) // '.stats.h5')
call note_file (trim(run%global%out_root) // '.resolved.toml')
call note_file (trim(run%global%out_root) // '.diag.txt')
call note_file (trim(run%global%out_root) // '.ledger.txt')
call note_file (trim(run%global%out_root) // '.import.txt')
call note_file (trim(run%global%out_root) // '.fp32.txt')
call note_file (trim(run%global%out_root) // '.migration.txt')
call fel_box_rule (fel_rule_bottom$, '', '', .false., .false.)

! What the result says of its own convergence, a line each. The power outside the
! split angle is the wide-angle emission of the point beamlets, and a slice is thin
! where it has under 128 beamlets or over a tenth of its power outside the mode.

any_warning = .false.
n_thin = 0
if (cvg%ok) then
  worst = maxval(cvg%ratio)
  if (worst > 0.1_rp) call warn_result ('outside / inside = ' // number_text(worst) // ', above 0.1')

  e_thin = 128
  if (run%bparam%load_mode == 'sample' .and. run%bparam%shot_noise .and. run%bparam%quiet_start) &
        e_thin = 128 * max(1.0_rp, run%bparam%beamlet_min_electrons)

  n_thin = 0;  n_shown = 0
  do is = 1, min(stats%nslice, size(run%fbeam%slice))

    ! A slice with no charge at the mode peak is not thin, and a slice with fewer
    ! electrons than e_thin has the beamlets its charge allows.

    if (stats%charge_live(is, cvg%ir(1)) <= 0) cycle
    nbl = run%fbeam%slice(is)%n / max(1, run%bparam%beamlet_size)
    ratio_s = 0
    if (stats%f_angles_valid(is, cvg%ir(1)) == 1 .and. stats%f_p_in(is, cvg%ir(1)) > 0) &
          ratio_s = max(0.0_rp, stats%f_power(is, cvg%ir(1)) - stats%f_p_in(is, cvg%ir(1))) / &
                    stats%f_p_in(is, cvg%ir(1))
    few = (nbl < 128 .and. stats%charge_live(is, cvg%ir(1)) >= e_thin * e_charge)
    if (.not. few .and. ratio_s <= 0.1_rp) cycle
    n_thin = n_thin + 1
    if (n_shown >= n_thin_max$) cycle
    n_shown = n_shown + 1
    write (line, '(a, i0, 3a, i0, 2a)') 'slice ', is, ' at ', &
          si_text((is - 1) * run%fbeam%slice_spacing, 'm'), ': ', nbl, ' beamlets, outside / inside ', &
          number_text(ratio_s)
    call warn_result (trim(line))
  enddo
  if (n_thin > n_shown) then
    write (line, '(a, i0, a, i0, a, i0, a)') 'and ', n_thin - n_shown, ' more: ', n_thin, ' thin slices of ', &
          stats%nslice, ' at the mode peak'
    call warn_result (trim(line))
  endif
endif

! The values of this printout and the list of files, into the stats file (fel_record_mod).

call fel_record_int ('run', 'n_thin_slice', 'The slices that are thin at the mode peak: fewer than 128 ' // &
      'beamlets, or more than 0.1 of the power outside the split angle.', n_thin)
call fel_record_write (trim(run%global%out_root) // '.stats.h5', .true., rec_err)

contains

!------------------------------------------------------------------------------
!+
! Function secs_right (seconds) result (txt)
!
! Routine to write seconds as fel_cost_secs_str does, set to the right of nine columns
! so that a column of times lines up.
!-

function secs_right (seconds) result (txt)

real(rp) seconds
character(9) txt
character(32) str

!

str = fel_cost_secs_str(seconds)
txt = adjustr(str(1:9))

end function secs_right

!------------------------------------------------------------------------------
!+
! Subroutine add_time (name, seconds)
!
! Routine to add one phase to the timing part, in the left column and the right in turn:
! its seconds and its share of the walk.
!-

subroutine add_time (name, seconds)

character(*) name
real(rp) seconds
character(40) val

!

write (val, '(a9, i5, a)') secs_right(seconds), nint(100 * seconds / walk), ' %'
n_left = n_left + 1
if (mod(n_left, 2) == 1) then
  call fel_column_add (c_left, name, trim(val))
else
  call fel_column_add (c_right, name, trim(val))
endif

end subroutine add_time

!------------------------------------------------------------------------------
!+
! Subroutine warn_result (text)
!
! Routine to write one warning line, with a blank line before the first.
!-

subroutine warn_result (text)

character(*) text

!

if (.not. any_warning) call out_io (s_blank$, r_name, '')
any_warning = .true.
call fel_print_line (' ! ' // text)

end subroutine warn_result


!+
! Subroutine note_file (name)
!
! Routine to list one output file with its size, if it is there. Existence is the test,
! so the footer never claims a file the run did not write.
!-

subroutine note_file (name)

character(*) name
integer(8) fsize
logical there

!

inquire (file = name, exist = there, size = fsize)
if (.not. there) return
call note_row (name, real(fsize, rp))

end subroutine note_file

!+
! Subroutine note_dumps (prefix, beam, field)
!
! Routine to list the beam dump <prefix>.beam.h5 and the field set's dumps
! <prefix>[-h<h>].wf.h5, each where its switch is on and the file is there.
!-

subroutine note_dumps (prefix, beam, field)

character(*) prefix
logical beam, field
integer jh
character(8) jsuf

!

if (beam) call note_file (prefix // '.beam.h5')
if (.not. field) return
do jh = 1, run%n_harm
  jsuf = ''
  if (run%ffield(jh)%harm /= 1) write (jsuf, '(a, i0)') '-h', run%ffield(jh)%harm
  call note_file (prefix // trim(jsuf) // '.wf.h5')
enddo

end subroutine note_dumps

!+
! Subroutine note_series (suffix)
!
! Routine to list the frame files <out_root>-<record><suffix> of records 1 to nrec as one
! row, the records as a shell brace range, with the size of the series. A series with
! one record is listed by its own name. Every file of the series must be there, since
! the frame writer writes one at every comb record, and a series with a file missing is
! not listed.
!-

subroutine note_series (suffix)

character(*) suffix
character(440) fname, pattern
integer(8) fsize
real(rp) total
integer jr
logical there

!

total = 0
do jr = 1, nrec
  write (fname, '(2a, i6.6, a)') trim(root), '-', jr, suffix
  inquire (file = trim(fname), exist = there, size = fsize)
  if (.not. there) return
  total = total + fsize
enddo
if (nrec == 1) then
  write (pattern, '(2a, i6.6, a)') trim(root), '-', 1, suffix
else
  write (pattern, '(2a, i6.6, a, i6.6, 2a)') trim(root), '-{', 1, '..', nrec, '}', suffix
endif
call note_row (trim(pattern), total)

end subroutine note_series

!+
! Subroutine note_row (name, fsize)
!
! Routine to write one row of the file list, the first after "Wrote".
!-

subroutine note_row (name, fsize)

character(*) name
character(500) row
character(12) size_text
real(rp) fsize

!

n_listed = n_listed + 1
call fel_record_file (name, fsize)
size_text = byte_text(fsize)
if (len_trim(name) <= 62) then
  row = name
  row(64:75) = adjustr(size_text)
else
  row = trim(name) // '   ' // trim(size_text)
endif
call fel_box_full (trim(row))

end subroutine note_row

!+
! Function record_name (text) result (name)
!
! Routine to make a dataset's name from a phase's printed name, the spaces as underscores.
!-

function record_name (text) result (name)

character(*) text
character(len(text)) name
integer i

!

name = text
do i = 1, len_trim(name)
  if (name(i:i) == ' ') name(i:i) = '_'
enddo

end function record_name

end subroutine fel_write_footer

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_read_openpmd_into_field (run, fname, ihh, err_flag)
!
! Routine to read an openPMD wavefront into field-set entry ihh (the fundamental
! import path). The photon energy must be the fundamental's. A file carrying a
! harmonic in field.file(1) is refused.
!
! Input:
!   run       -- fel_run_struct: Run state.
!   fname     -- character(*): openPMD wavefront file to read.
!   ihh       -- integer: Field-set index to fill.
!
! Output:
!   run       -- fel_run_struct: run%ffield(ihh)%wf holds the imported field.
!   err_flag  -- logical: Set True on a read error or wavelength mismatch. False otherwise.
!-

subroutine fel_read_openpmd_into_field (run, fname, ihh, err_flag)

type (fel_run_struct), target :: run
type (fel_field_struct), pointer :: ffield(:)
type (fel_beam_struct), pointer :: fbeam
type (hdf5_info_struct) info
integer(hid_t) f_id, g_id
character(*) fname
integer ihh, ihc, h5_err
real(rp) e_photon, residual, s_c
integer n_x, n_y
logical rerr, err_flag, has_res, dark
character(60) fmt_str
character(40) name_c
character(*), parameter :: r_name = 'fel_read_openpmd_into_field'

!

err_flag = .true.
ffield => run%ffield
fbeam => run%fbeam

call wavefront_read_openpmd (ffield(ihh)%wf, fname, rerr, e_photon)
if (rerr) return
if (abs(ffield(ihh)%wf%wavelength - fbeam%wavelength) > 1e-6_rp * fbeam%wavelength) then
  call out_io (s_error$, r_name, 'THE openPMD FILE IN FIELD_FILE(1) DOES NOT CARRY THE FUNDAMENTAL:', &
               'ITS PHOTONENERGY WAVELENGTH VS THE BEAM: \2es20.12\ ', &
               r_array = [ffield(ihh)%wf%wavelength, fbeam%wavelength])
  err_flag = .true.;  return
endif
! The file's planes are the field's own spacing apart, the beam's slice times the
! slices a plane spans (fel_field_sample), to the rounding a writer leaves. A dark file
! at another spacing carries nothing to resample, so its transverse grid alone is taken
! and the record is built dark on the field's own planes: a SASE input file whose field
! file fixes the grid then runs at any sampling.

if (abs(ffield(ihh)%wf%dz - ffield(ihh)%n_per * fbeam%slice_spacing) > 1e-9_rp * fbeam%slice_spacing) then
  dark = all(ffield(ihh)%wf%Ex == 0)
  if (dark .and. allocated(ffield(ihh)%wf%Ey)) dark = all(ffield(ihh)%wf%Ey == 0)
  if (.not. dark) then
    call out_io (s_error$, r_name, 'THE FIELD FILE''S PLANES ARE \es14.6\ m APART, AND THE RUN''S FIELD', &
                 'SAMPLES EVERY \es14.6\ m: \i0\ WAVELENGTHS, \i0\ BEAM SLICES A PLANE. THE FILE: ' // trim(fname), &
                 r_array = [ffield(ihh)%wf%dz, ffield(ihh)%n_per * fbeam%slice_spacing], &
                 i_array = [ffield(ihh)%n_wavelength, ffield(ihh)%n_per])
    return
  endif
  n_x = size(ffield(ihh)%wf%Ex, 1);  n_y = size(ffield(ihh)%wf%Ex, 2)
  deallocate (ffield(ihh)%wf%Ex)
  allocate (ffield(ihh)%wf%Ex(n_x, n_y, run%nfield))
  ffield(ihh)%wf%Ex = 0
  if (allocated(ffield(ihh)%wf%Ey)) then
    deallocate (ffield(ihh)%wf%Ey)
    allocate (ffield(ihh)%wf%Ey(n_x, n_y, run%nfield))
    ffield(ihh)%wf%Ey = 0
  endif
  call out_io (s_info$, r_name, 'The field file is dark and its planes are \es12.4\ m apart, so its grid alone is ' // &
               'taken and the field starts dark on \i0\ planes \es12.4\ m apart.', &
               r_array = [ffield(ihh)%wf%dz, ffield(ihh)%n_per * fbeam%slice_spacing], i_array = [run%nfield])
endif
ffield(ihh)%wf%dz = ffield(ihh)%n_per * fbeam%slice_spacing
ffield(ihh)%wf%wavelength = fbeam%wavelength   ! One wavelength for beam and field, as the
                                               !   1e-12 consistency check requires.

! Read the slippage the record had accumulated and not yet rotated when it was written
! (fel_write_slippage_residual). The record itself does not hold it. Also read the
! checkpoint group saying where the frame was taken. A continuation requires both. It
! also requires the position and element to match the beam file's, since a beam and a
! field from two boundaries are two states and not one.
!
! An initialization takes neither. The record starts from a whole-slice boundary. A
! file that carried a residual gets a message saying so, so a stale residual never
! enters a fresh run. A residual the tracker could not have written is refused in either
! mode. Its magnitude never reaches n_wavelength, since the threshold is 0.8 of it.

call hdf5_open_file (fname, 'READ', f_id, rerr, .false.)
if (rerr) return
g_id = fel_private_group (f_id, 'field', .false., rerr)
if (rerr) then
  call h5fclose_f (f_id, h5_err)
  return
endif
has_res = .false.
residual = 0
if (g_id >= 0) then
  info = hdf5_attribute_info (g_id, 'slippageResidual', rerr, .false.)
  has_res = .not. rerr
  if (has_res) call hdf5_read_attribute_real (g_id, 'slippageResidual', residual, rerr, .true.)
  call H5Gclose_f (g_id, h5_err)
  if (has_res .and. rerr) then
    call h5fclose_f (f_id, h5_err)
    return
  endif
endif
if (has_res) then
  if (.not. ieee_is_finite(residual) .or. abs(residual) >= ffield(ihh)%n_wavelength) then
    call h5fclose_f (f_id, h5_err)
    call out_io (s_error$, r_name, 'THE FIELD FILE''S slippageResidual IS NOT A VALUE THIS TRACKER WRITES: \es12.4\ ', &
                 'ITS MAGNITUDE MUST STAY BELOW THE FIELD''S n_wavelength = \i0\. THE FILE: ' // trim(fname), &
                 r_array = [residual], i_array = [ffield(ihh)%n_wavelength])
    return
  endif
endif

if (.not. run%global%continuation) then
  call h5fclose_f (f_id, h5_err)
  if (has_res) then
    call out_io (s_info$, r_name, 'The field file carries a slippage residual of \es12.4\ wavelengths, ' // &
                 'ignored: run.continuation is off, so the record starts from a whole-slice boundary.', &
                 r_array = [residual])
  else
    call out_io (s_info$, r_name, 'The field file carries no slippage residual, so the record starts from a ' // &
                 'whole-slice boundary.')
  endif
  err_flag = .false.
  return
endif

if (.not. has_res) then
  call h5fclose_f (f_id, h5_err)
  call out_io (s_error$, r_name, 'CONTINUATION FROM A FIELD FILE WITH NO slippageResidual: ' // trim(fname), &
               'THE RESIDUAL IS PART OF THE FIELD''S STATE, AND A CONTINUATION TAKES NOTHING FROM A FILE THAT LACKS IT.')
  return
endif

g_id = fel_private_group (f_id, 'checkpoint', .false., rerr)
if (rerr) then
  call h5fclose_f (f_id, h5_err)
  return
endif
if (g_id < 0) then
  call h5fclose_f (f_id, h5_err)
  call out_io (s_error$, r_name, 'CONTINUATION FROM A FIELD FILE WITH NO CHECKPOINT GROUP: ' // trim(fname), &
       'A CHECKPOINT IS WRITTEN AT AN ELEMENT BOUNDARY OUTSIDE EVERY UNDULATOR OR AT THE PHYSICAL END', &
       'OF ONE. A FRAME FROM INSIDE AN ELEMENT, FROM A SLICED ENDPOINT INSIDE AN UNDULATOR, OR FROM', &
       'BEFORE THE GROUP EXISTED IS NOT ONE, AND NOTHING IS TAKEN FROM IT.')
  return
endif
fmt_str = ''          ! The reads fill the characters they have and leave the rest.
name_c = ''
call hdf5_read_attribute_string (g_id, 'checkpointFormat', fmt_str, rerr, .true.)
if (.not. rerr) call hdf5_read_attribute_real (g_id, 'sPosition', s_c, rerr, .true.)
if (.not. rerr) call hdf5_read_attribute_string (g_id, 'elementName', name_c, rerr, .true.)
call H5Gclose_f (g_id, h5_err)
call h5fclose_f (f_id, h5_err)
if (rerr) then
  call out_io (s_error$, r_name, 'CHECKPOINT MEMBER MISSING FROM THE FIELD FILE''S GROUP: ' // trim(fname), &
               'THE GROUP IS INCOMPLETE AND NOTHING IS TAKEN FROM IT.')
  return
endif
if (trim(fmt_str) /= fel_checkpoint_format$) then
  do ihc = 1, len(fmt_str)          ! Bytes that are not text print as question marks.
    if (iachar(fmt_str(ihc:ihc)) < 32 .or. iachar(fmt_str(ihc:ihc)) > 126) fmt_str(ihc:ihc) = '?'
  enddo
  call out_io (s_error$, r_name, 'CHECKPOINT VERSION NOT KNOWN TO THIS TRACKER: "' // trim(fmt_str) // '".', &
       'THIS TRACKER READS "' // fel_checkpoint_format$ // '". THE FILE: ' // trim(fname))
  return
endif
if (.not. run%ckpt%loaded) then
  call out_io (s_error$, r_name, 'THE FIELD FILE WAS READ BEFORE THE BEAM FILE''S CHECKPOINT. PLEASE REPORT THIS!')
  return
endif
if (abs(s_c - run%ckpt%s) > 1.0e-9_rp .or. trim(name_c) /= trim(run%ckpt%ele_name)) then
  call out_io (s_error$, r_name, 'CHECKPOINT POSITIONS DISAGREE: THE BEAM FILE WAS WRITTEN AT s = \es16.8\ m ' // &
       'AFTER ELEMENT ' // trim(run%ckpt%ele_name) // ',', &
       'AND THE FIELD FILE AT s = \es16.8\ m AFTER ELEMENT ' // trim(name_c) // '. A CONTINUATION TAKES', &
       'ONE BOUNDARY''S PAIR OF FRAMES. THE FIELD FILE: ' // trim(fname), r_array = [run%ckpt%s, s_c])
  return
endif

ffield(ihh)%slip%accuslip = residual
call out_io (s_info$, r_name, 'The field file carries a slippage residual of \es12.4\ wavelengths, restored, ' // &
             'and its checkpoint agrees with the beam file''s.', r_array = [residual])

err_flag = .false.

end subroutine fel_read_openpmd_into_field

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_import_harmonic_field (run, fname, err_flag)
!
! Routine to import a harmonic field. The file's photon energy is matched to the
! field-set entry carrying that harmonic, and no match is refused. The grid and window
! must be the fundamental's. The walk's wavelength convention (fundamental / h) is kept.
!
! Input:
!   run       -- fel_run_struct: Run state with the fundamental already loaded.
!   fname     -- character(*): openPMD wavefront file to read.
!
! Output:
!   run       -- fel_run_struct: The matching run%ffield entry holds the imported field.
!   err_flag  -- logical: Set True on a read error, no matching harmonic, or a grid
!                  mismatch. False otherwise.
!-

subroutine fel_import_harmonic_field (run, fname, err_flag)

type (fel_run_struct), target :: run
type (fel_field_struct), pointer :: ffield(:)
type (fel_beam_struct), pointer :: fbeam
type (wavefront_struct), pointer :: wf
character(*) fname
type (wavefront_struct) wtmp
real(rp) e_photon, e1
integer ihh, n_harm
logical rerr, err_flag
character(*), parameter :: r_name = 'fel_import_harmonic_field'

!

err_flag = .true.
ffield => run%ffield
fbeam => run%fbeam
wf => run%ffield(1)%wf
n_harm = run%n_harm

call wavefront_read_openpmd (wtmp, fname, rerr, e_photon)
if (rerr) return
e1 = h_planck * c_light / fbeam%wavelength * e_charge

do ihh = 2, n_harm
  if (abs(e_photon - ffield(ihh)%harm * e1) <= 1e-6_rp * ffield(ihh)%harm * e1) exit
enddo
if (ihh > n_harm) then
  call out_io (s_error$, r_name, 'THE PHOTONENERGY OF ' // trim(fname), &
               '(\es13.5\ J) MATCHES NO FIELD OF THIS RUN''S HARMONICS LIST.', r_array = [e_photon])
  err_flag = .true.;  return
endif

if (size(wtmp%Ex,1) /= size(wf%Ex,1) .or. size(wtmp%Ex,3) /= size(wf%Ex,3) .or. &
    abs(wtmp%dx - wf%dx) > 1e-12_rp * wf%dx) then
  call out_io (s_error$, r_name, 'HARMONIC IMPORT ' // trim(fname), &
               'DOES NOT MATCH THE FUNDAMENTAL''S GRID AND WINDOW (ONE TIME WINDOW, ONE GRID).')
  err_flag = .true.;  return
endif

call move_alloc (wtmp%Ex, ffield(ihh)%wf%Ex)
if (allocated(wtmp%Ey)) call move_alloc (wtmp%Ey, ffield(ihh)%wf%Ey)
ffield(ihh)%wf%dx = wtmp%dx;  ffield(ihh)%wf%dy = wtmp%dy
ffield(ihh)%wf%dz = fbeam%slice_spacing
ffield(ihh)%wf%wavelength = wf%wavelength / ffield(ihh)%harm

err_flag = .false.

end subroutine fel_import_harmonic_field

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_drain_bank (run, ihh, err_flag)
!
! Routine to drain the slippage bank. Each slice carried out streams to the pulse file,
! <out_root>-pulse.wf.h5, as it leaves, so peak memory is a handful of grid planes and
! the pulse is the one file on disk. With write_escaped it also streams to the escaped
! file, <out_root>-escaped.wf.h5. Both are openPMD EXT_Wavefront files
! (wavefront_openpmd_mod), made with the live window's grid and wavelength and
! component datasets that grow along z (wavefront_openpmd_create). Their planes are in
! V/m, as the dumps are.
!
! The files are made when the first plane with light leaves. The planes before it are
! the empty slices ahead of the bunch, all zero, and they are written then as the
! files' first banked planes, so the files are the same as if they had been made at the
! first slice. A run whose light never leaves the window, as a derived window holds the
! slippage, makes no file at all. Until the files are made, run%n_banked counts the
! zero planes.
!
! Each slice's wavefront_params and transmission z are kept in run%bank_pms and
! run%bank_z, the params the full 4x4 at bank time, which the analytic free-space
! propagation of pulse statistics needs. They go into the files at finalize
! (fel_write_bank_group).
!
! During the run the banked planes are in transmission order: banked plane k is the
! k-th to leave. That is the reverse of the mesh's order. Light leaves at the head, so
! the first to leave is the furthest ahead. In the pulse they follow nslice planes the
! window fills at finalize, and fel_finish_pulse reverses them there as it drifts them
! to the exit. In the escaped file they are its only planes, and fel_order_bank
! reverses them.
!
! No one beamline position holds the escaped planes, since each is as it was when it
! left at its own z. So the escaped file's zCoordinate is NaN, and
! lucifer/escaped/zTransmit carries each plane's own. Its mesh's z origin is the slot
! just past the window's head, nslice slice spacings along the pulse. The pulse's
! zCoordinate is written at finalize, as the exit plane.
!
! Input:
!   run       -- fel_run_struct: Run state. run%ffield(ihh)%bank holds the slices to drain.
!   ihh       -- integer: Field-set index whose bank is drained.
!
! Output:
!   run       -- fel_run_struct: Banked params appended to run%bank_z, run%bank_pms and,
!                  with a live y component, run%bank_pms_y. run%n_banked(ihh),
!                  run%pul_id(ihh) and run%esc_id(ihh) updated.
!   err_flag  -- logical: Set True on a write error. False otherwise.
!-

subroutine fel_drain_bank (run, ihh, err_flag)

type (fel_run_struct), target :: run
type (fel_field_struct), pointer :: ffield(:)
type (fel_beam_struct), pointer :: fbeam
real(rp), pointer :: z_now
logical err_flag
integer n_harm

type (wavefront_params_struct) pms
type (fel_bank_struct), pointer :: bnk
type (wavefront_struct), pointer :: wfl
integer(hid_t) m_id
real(rp), allocatable :: gz(:,:), gp(:,:,:), gpy(:,:,:)
complex(wf_rp), allocatable :: zero(:,:,:)
integer ihh, k, n0, n_cap, h5e
logical berr, two, dark, fresh

!

err_flag = .true.
ffield => run%ffield
fbeam => run%fbeam
z_now => run%z_now
n_harm = run%n_harm

bnk => ffield(ihh)%bank
wfl => ffield(ihh)%wf
if (bnk%n == 0) then
  err_flag = .false.
  return
endif
two = allocated(wfl%Ey)
n0 = run%n_banked(ihh)

! The params, kept in memory in transmission order. The arrays double when full.

if (.not. allocated(run%bank_z)) then
  n_cap = max(1024, n0 + bnk%n)
  allocate (run%bank_z(n_cap, n_harm), run%bank_pms(25, n_cap, n_harm))
  if (two) allocate (run%bank_pms_y(25, n_cap, n_harm))
elseif (n0 + bnk%n > size(run%bank_z, 1)) then
  n_cap = max(2 * size(run%bank_z, 1), n0 + bnk%n)
  call move_alloc (run%bank_z, gz);  call move_alloc (run%bank_pms, gp)
  allocate (run%bank_z(n_cap, n_harm), run%bank_pms(25, n_cap, n_harm))
  run%bank_z(1:size(gz,1), :) = gz;  run%bank_pms(:, 1:size(gz,1), :) = gp
  if (two) then
    call move_alloc (run%bank_pms_y, gpy)
    allocate (run%bank_pms_y(25, n_cap, n_harm))
    run%bank_pms_y(:, 1:size(gz,1), :) = gpy
  endif
endif

do k = 1, bnk%n
  call wavefront_params_of_plane (bnk%plane(:,:,k), wfl%dx, wfl%wavelength, &
                                  wfl%dz, pms, .true., berr)
  if (berr) return
  pms%s = z_now
  run%bank_z(n0 + k, ihh) = z_now
  run%bank_pms(:, n0 + k, ihh) = [pms%centroid, reshape(pms%sigma, [16]), &
                           pms%energy, pms%power, pms%on_axis_intensity, pms%emit_x, pms%emit_y]
  if (two) then
    call wavefront_params_of_plane (bnk%plane_y(:,:,k), wfl%dx, wfl%wavelength, &
                                    wfl%dz, pms, .true., berr)
    if (berr) return
    run%bank_pms_y(:, n0 + k, ihh) = [pms%centroid, reshape(pms%sigma, [16]), &
                           pms%energy, pms%power, pms%on_axis_intensity, pms%emit_x, pms%emit_y]
  endif
enddo

! The files, at the first plane with light. The pulse's first nslice planes are the
! window's, written at finalize.

fresh = .false.
if (run%pul_id(ihh) == 0) then
  dark = all(bnk%plane(:,:,1:bnk%n) == 0)
  if (dark .and. two) dark = all(bnk%plane_y(:,:,1:bnk%n) == 0)
  if (dark) then
    run%n_banked(ihh) = n0 + bnk%n
    err_flag = .false.
    return
  endif
  call wavefront_openpmd_create (trim(fel_bank_file_name(run, ihh, 'pulse')), wfl, run%nfield, &
                                 run%pul_id(ihh), m_id, berr, extendable = .true.)
  if (berr) return
  call H5Gclose_f (m_id, h5e)
  if (run%global%write_escaped) then
    call wavefront_openpmd_create (trim(fel_bank_file_name(run, ihh, 'escaped')), wfl, 0, run%esc_id(ihh), m_id, &
                                   berr, z_offset = run%nfield * wfl%dz, extendable = .true., &
                                   z_coordinate = ieee_value(0.0_rp, ieee_quiet_nan))
    if (berr) return
    call H5Gclose_f (m_id, h5e)
  endif
  fresh = .true.
endif

call append_planes (run%pul_id(ihh), run%nfield);  if (berr) return
if (run%esc_id(ihh) /= 0) call append_planes (run%esc_id(ihh), 0)
if (berr) return
run%n_banked(ihh) = n0 + bnk%n

err_flag = .false.

!------------------------------------------------------------------------------
contains

!+
! Subroutine append_planes (f_id, i_first)
!
! Routine to write this drain's planes into one open file after the banked planes it
! holds, which follow i_first planes, and in a new file after the zero planes that left
! before it. Sets berr on a write error.
!-

subroutine append_planes (f_id, i_first)

integer(hid_t) f_id
integer i_first

!

call wavefront_openpmd_mesh (f_id, m_id, berr);  if (berr) return
if (fresh .and. n0 > 0) then
  if (.not. allocated(zero)) allocate (zero(size(wfl%Ex, 1), size(wfl%Ex, 2), 1))
  zero = 0
  do k = 1, n0
    call wavefront_openpmd_write_planes (m_id, 'x', i_first + k, zero, berr)
    if (.not. berr .and. two) call wavefront_openpmd_write_planes (m_id, 'y', i_first + k, zero, berr)
    if (berr) exit
  enddo
endif
if (.not. berr) call wavefront_openpmd_write_planes (m_id, 'x', i_first + n0 + 1, bnk%plane(:,:,1:bnk%n), berr)
if (.not. berr .and. two) call wavefront_openpmd_write_planes (m_id, 'y', i_first + n0 + 1, &
                                                               bnk%plane_y(:,:,1:bnk%n), berr)
call H5Gclose_f (m_id, h5e)

end subroutine append_planes

end subroutine fel_drain_bank

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_order_bank (run, ihh, err_flag)
!
! Routine to put the escaped file in order at finalize, once n_banked is known.
! fel_drain_bank writes the planes in transmission order as they leave. Light leaves at
! the window's head, so the k-th to leave is the (n_banked - k + 1)-th in increasing
! position along the pulse, which is the mesh's order. The planes are reversed in place,
! two at a time, so peak memory stays at a few planes. The group lucifer/escaped is then
! written (fel_write_bank_group). run%bank_z and run%bank_pms stay in transmission order.
!
! Only a run with write_escaped has an escaped file. The pulse is built without it
! (fel_finish_pulse).
!
! Input:
!   run       -- fel_run_struct: Run state. run%esc_id(ihh) is the open escaped file.
!   ihh       -- integer: Field-set index.
!
! Output:
!   err_flag  -- logical: Set True on an HDF5 error. False otherwise.
!-

subroutine fel_order_bank (run, ihh, err_flag)

type (fel_run_struct), target :: run
type (wavefront_struct), pointer :: wfl
complex(wf_rp), allocatable :: lo(:,:,:), hi(:,:,:)
integer(hid_t) m_id
integer ihh, k, n, h5e
logical err_flag, berr, two
character(*), parameter :: r_name = 'fel_order_bank'

!

err_flag = .true.
n = run%n_banked(ihh)
wfl => run%ffield(ihh)%wf
two = allocated(wfl%Ey)

call wavefront_openpmd_mesh (run%esc_id(ihh), m_id, berr);  if (berr) return
allocate (lo(size(wfl%Ex, 1), size(wfl%Ex, 2), 1), hi(size(wfl%Ex, 1), size(wfl%Ex, 2), 1))
berr = .false.
do k = 1, n / 2
  call swap_planes ('x');  if (berr) exit
  if (two) call swap_planes ('y')
  if (berr) exit
enddo
call H5Gclose_f (m_id, h5e)
if (berr) then
  call out_io (s_error$, r_name, 'CANNOT REORDER THE ESCAPED FILE ' // trim(fel_bank_file_name(run, ihh, 'escaped')))
  return
endif

call fel_write_bank_group (run%esc_id(ihh), run, ihh, 'escaped', &
      'Each plane of the mesh left the window at its own beamline position, zTransmit, and is ' // &
      'as it was then. No one position holds them all, so zCoordinate is NaN.', err_flag)

!------------------------------------------------------------------------------
contains

!+
! Subroutine swap_planes (name)
!
! Routine to exchange planes k and n - k + 1 of one component.
!-

subroutine swap_planes (name)

character(*) name

!

call wavefront_openpmd_read_planes (m_id, name, k, lo, berr);  if (berr) return
call wavefront_openpmd_read_planes (m_id, name, n - k + 1, hi, berr);  if (berr) return
call wavefront_openpmd_write_planes (m_id, name, k, hi, berr);  if (berr) return
call wavefront_openpmd_write_planes (m_id, name, n - k + 1, lo, berr)

end subroutine swap_planes

end subroutine fel_order_bank

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_write_bank_group (f_id, run, ihh, group, description, err_flag)
!
! Routine to write what openPMD has no name for beside a file's standard records: where
! each banked plane left the window, and its moments then. The group is
! lucifer/<group>, lucifer/escaped in the escaped file and lucifer/pulse in the pulse,
! and holds zTransmit and wavefrontParams, with wavefrontParamsY for a y component, one
! row per banked plane in the mesh's order. That is the reverse of run%bank_z's
! transmission order. The two files' groups hold the same rows.
!
! Input:
!   f_id        -- integer(hid_t): The open file.
!   run         -- fel_run_struct: Run state. run%bank_z and run%bank_pms hold the rows.
!   ihh         -- integer: Field-set index.
!   group       -- character(*): 'escaped' or 'pulse'.
!   description -- character(*): The group's description attribute.
!
! Output:
!   err_flag    -- logical: Set True on an HDF5 error. False otherwise.
!-

subroutine fel_write_bank_group (f_id, run, ihh, group, description, err_flag)

type (fel_run_struct), target :: run
integer(hid_t) f_id, g_id
integer ihh, n, h5e
logical err_flag, berr, bad, two
character(*) group, description
! The 25 columns of wavefrontParams: wavefront_params_struct flattened, its sigma(4,4)
! column by column, over the axes (x, theta_x, y, theta_y).
character(20), parameter :: columns(25) = [character(20) :: &
      'centroid_x', 'centroid_theta_x', 'centroid_y', 'centroid_theta_y', &
      'sigma_11', 'sigma_21', 'sigma_31', 'sigma_41', 'sigma_12', 'sigma_22', 'sigma_32', 'sigma_42', &
      'sigma_13', 'sigma_23', 'sigma_33', 'sigma_43', 'sigma_14', 'sigma_24', 'sigma_34', 'sigma_44', &
      'energy', 'power', 'on_axis_intensity', 'emit_x', 'emit_y']
character(*), parameter :: params_descrip = &
      'wavefront_params_struct of each slice as it left the window, in the mesh''s order: ' // &
      'the intensity-weighted centroid (x, theta_x, y, theta_y) [m, rad], the Wigner second-moment ' // &
      'matrix over those axes column by column [m^2, m rad, rad^2], the energy [J], the power [W], ' // &
      'the on-axis intensity [W/m^2] and the emittances [m rad]. The columns attribute names them.'

!

err_flag = .true.
n = run%n_banked(ihh)
two = allocated(run%ffield(ihh)%wf%Ey)

g_id = fel_private_group (f_id, group, .true., berr);  if (berr) return
bad = .false.
call hdf5_write_attribute_string (g_id, 'description', description, berr)
bad = bad .or. berr
call fel_h5_real (g_id, 'zTransmit', 'm', 'transmission position', &
      'The beamline position at which each banked slice left the window, in the mesh''s order.', &
      '', run%bank_z(n:1:-1, ihh), bad)
call fel_h5_real (g_id, 'wavefrontParams', 'mixed', 'wavefront params', params_descrip, &
      '', run%bank_pms(:, n:1:-1, ihh), bad)
call fel_h5_dset_attr_strs (g_id, 'wavefrontParams', 'columns', columns, bad)
if (two) then
  call fel_h5_real (g_id, 'wavefrontParamsY', 'mixed', 'wavefront params, y', &
        'The y component''s wavefrontParams, the same columns.', '', run%bank_pms_y(:, n:1:-1, ihh), bad)
  call fel_h5_dset_attr_strs (g_id, 'wavefrontParamsY', 'columns', columns, bad)
endif
call H5Gclose_f (g_id, h5e)
if (bad) return

err_flag = .false.

end subroutine fel_write_bank_group

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_bank_file_name (run, ihh, what) result (fname)
!
! Routine to construct the escaped or the pulse file name of field ihh,
! <out_root>-escaped.wf.h5 or <out_root>-pulse.wf.h5. The fundamental keeps its
! pre-harmonic name. A harmonic's file carries -h<h>.
!
! Input:
!   run    -- fel_run_struct: Run state (out_root and the field set).
!   ihh    -- integer: Field-set index.
!   what   -- character(*): 'escaped' or 'pulse'.
!
! Output:
!   fname  -- character(420): The file name.
!-

function fel_bank_file_name (run, ihh, what) result (fname)

type (fel_run_struct), target :: run
type (fel_field_struct), pointer :: ffield(:)
character(*) what
character(400) out_root
integer ihh
character(420) fname

ffield => run%ffield
out_root = run%global%out_root

if (ffield(ihh)%harm == 1) then
  fname = trim(out_root) // '-' // what // '.wf.h5'
else
  write (fname, '(4a, i0, a)') trim(out_root), '-', what, '-h', ffield(ihh)%harm, '.wf.h5'
endif

end function fel_bank_file_name

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_write_stats_file (run, err_flag)
!
! Routine to write the stats file <out_root>.stats.h5: the records taken, the input as
! resolved in params/, the lattice table, meta/ and the values of the printout's header
! in setup/ (fel_record_mod). A load_only run calls it with no record taken, and its
! file then holds what setup worked out and no tracking.
!
! Input:
!   run      -- fel_run_struct: Run state.
!
! Output:
!   err_flag -- logical: Set true if the file could not be written.
!-

subroutine fel_write_stats_file (run, err_flag)

type (fel_run_struct), target :: run
type (fel_stats_params_struct) sprm
character(:), allocatable :: file_name
logical err_flag

! The parameters the file states as data. They are assembled here because the whole run
! is visible here. fel_stats_mod knows the accumulator and nothing else.

file_name = trim(run%global%out_root) // '.stats.h5'
sprm%slice_spacing = run%fbeam%slice_spacing
sprm%species = 'electron'
sprm%beta0 = fel_p0_mc(run%fbeam) / sqrt(fel_p0_mc(run%fbeam)**2 + 1)
sprm%grid_dx = run%ffield(1)%wf%dx
sprm%grid_dy = run%ffield(1)%wf%dy

call fel_stats_write (run%stats, sprm, file_name, err_flag)
call fel_write_params (run, file_name)
if (err_flag) return
call fel_write_lattice (run, file_name)
call fel_write_meta (run, file_name)
call fel_record_write (file_name, .false., err_flag)

end subroutine fel_write_stats_file

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_finalize_diagnostics (run, err_flag)
!
! Routine to finalize the run's file diagnostics: the stats file always. With
! keep_escaped_field, also the full pulse at the exit plane (fel_finish_pulse), and with
! write_escaped the escaped file's order and private group (fel_order_bank).
!
! The pulse's banked planes streamed into it during the run (fel_drain_bank). Each is
! free-space propagated over z_end - z_transmit here. Transmitted light is fixed
! information, and undulator vacuum is free space for light with no beam under it. The
! banked planes follow the live window's. Earliest-transmitted light is furthest ahead.
! So banked slice k lands at pulse index nslice + (run%n_banked - k + 1), which is also
! its plane in a kept escaped file. The caller has already unrotated the live window
! (the final-dump cshift).
!
! Both files are openPMD EXT_Wavefront files whose planes run along z in increasing
! position, as the dumps' do. The z of gridGlobalOffset is the position along the pulse
! of a file's first plane: 0 for the pulse, whose first plane is the live window's tail,
! and nslice times the slice spacing for the escaped file, whose first plane is the
! slot just past the window's head. Neither is the beamline position, which is
! zCoordinate: the exit plane for the pulse, and NaN for the escaped file, whose planes
! each left at their own (lucifer/escaped/zTransmit). The pulse carries the same rows in
! lucifer/pulse, so the one file holds all the escaped file holds but the raw planes.
!
! Input:
!   run       -- fel_run_struct: Run state at the end of tracking.
!
! Output:
!   run       -- fel_run_struct: Pulse and escaped files closed (run%pul_id and
!                  run%esc_id zeroed).
!   err_flag  -- logical: Set True on a write error. False otherwise.
!-

subroutine fel_finalize_diagnostics (run, err_flag)

type (fel_run_struct), target :: run
type (fel_field_struct), pointer :: ffield(:)
type (fel_beam_struct), pointer :: fbeam
integer h5e, ihh, n_harm, iu_del, ios_del
logical ferr, err_flag, keep_escaped_field
character(400) out_root
character(*), parameter :: r_name = 'fel_finalize_diagnostics'

!

err_flag = .true.
ffield => run%ffield
fbeam => run%fbeam
n_harm = run%n_harm
keep_escaped_field = run%global%keep_escaped_field
out_root = run%global%out_root

call fel_write_stats_file (run, ferr)
if (ferr) return

if (.not. keep_escaped_field) then
  err_flag = .false.
  return
endif

do ihh = 1, n_harm

  ! A bank that holds no light, only the empty slices that left ahead of the bunch, is
  ! no output, since the pulse is then the window, which the final dump holds.
  ! fel_drain_bank makes no file for such a bank, and files that hold no energy are
  ! removed.
  if (run%pul_id(ihh) == 0) then
    run%n_banked(ihh) = 0
    cycle
  endif
  if (ffield(ihh)%slip%u_escaped == 0) then
    call H5Fclose_f (run%pul_id(ihh), h5e)
    run%pul_id(ihh) = 0
    call delete_file (fel_bank_file_name(run, ihh, 'pulse'))
    if (run%esc_id(ihh) /= 0) then
      call H5Fclose_f (run%esc_id(ihh), h5e)
      run%esc_id(ihh) = 0
      call delete_file (fel_bank_file_name(run, ihh, 'escaped'))
    endif
    run%n_banked(ihh) = 0
    cycle
  endif

  if (run%esc_id(ihh) /= 0) then
    call fel_order_bank (run, ihh, ferr)
    if (ferr) return
    call H5Fclose_f (run%esc_id(ihh), h5e)
    run%esc_id(ihh) = 0
    call out_io (s_info$, r_name, 'Wrote ' // trim(fel_bank_file_name(run, ihh, 'escaped')))
  endif

  call fel_finish_pulse (run, ihh, ferr)
  if (ferr) return
  call out_io (s_info$, r_name, 'Wrote ' // trim(fel_bank_file_name(run, ihh, 'pulse')))
enddo

err_flag = .false.

!------------------------------------------------------------------------------
contains

!+
! Subroutine delete_file (name)
!
! Routine to remove one file of the bank.
!-

subroutine delete_file (name)

character(*) name

!

open (newunit = iu_del, file = trim(name), status = 'old', iostat = ios_del)
if (ios_del == 0) close (iu_del, status = 'delete')

end subroutine delete_file

end subroutine fel_finalize_diagnostics

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_finish_pulse (run, ihh, err_flag)
!
! Routine to finish the pulse file of one field, which fel_drain_bank opened and filled
! with the banked planes in transmission order after nslice empty planes. Banked plane
! k left at run%bank_z(k) and is drifted to the exit over z_end - run%bank_z(k). Its
! place in the mesh's order is nslice + n_banked - k + 1, so planes k and
! n_banked - k + 1 are read, drifted and exchanged two at a time, and peak memory is a
! few planes whatever the pulse's length. The live window then fills planes 1 to nslice,
! the group lucifer/pulse takes each banked plane's transmission z and params
! (fel_write_bank_group), and zCoordinate is set to the exit plane, the window's
! position, which the walk set to the exit. The file is then closed.
!
! Input:
!   run       -- fel_run_struct: Run state. run%pul_id(ihh) is the open pulse file.
!   ihh       -- integer: Field-set index.
!
! Output:
!   run       -- fel_run_struct: run%pul_id(ihh) closed and zeroed.
!   err_flag  -- logical: Set True on a read or write error. False otherwise.
!-

subroutine fel_finish_pulse (run, ihh, err_flag)

type (fel_run_struct), target :: run
type (wavefront_struct) lo, hi
type (wavefront_struct), pointer :: wfl
logical err_flag
integer(hid_t) pm_id
integer ihh, k, n, nslice, h5e
logical ferr, two
character(300) descrip
character(*), parameter :: r_name = 'fel_finish_pulse'

!

err_flag = .true.
nslice = run%nfield       ! The pulse's window is the field's planes.
n = run%n_banked(ihh)
wfl => run%ffield(ihh)%wf
two = allocated(wfl%Ey)

! Two one-plane wavefronts at the window's grid and wavelength. Ex and Ey do not couple
! in free space, so the one drift carries both.

call one_plane (lo);  call one_plane (hi)

call wavefront_openpmd_mesh (run%pul_id(ihh), pm_id, ferr);  if (ferr) return
do k = 1, (n + 1) / 2
  call read_plane (lo, k);  if (ferr) exit
  call wavefront_drift (lo, run%z_now - run%bank_z(k, ihh), ferr);  if (ferr) exit
  if (k /= n - k + 1) then
    call read_plane (hi, n - k + 1);  if (ferr) exit
    call wavefront_drift (hi, run%z_now - run%bank_z(n - k + 1, ihh), ferr);  if (ferr) exit
    call write_plane (hi, k);  if (ferr) exit
  endif
  call write_plane (lo, n - k + 1);  if (ferr) exit
enddo

! The live window ahead of them, and the exit plane as the file's beamline position.

if (.not. ferr) call wavefront_openpmd_write_planes (pm_id, 'x', 1, wfl%Ex, ferr)
if (.not. ferr .and. two) call wavefront_openpmd_write_planes (pm_id, 'y', 1, wfl%Ey, ferr)
if (.not. ferr) call hdf5_write_attribute_real (pm_id, 'zCoordinate', wfl%s_position, ferr)
call H5Gclose_f (pm_id, h5e)
if (ferr) then
  call out_io (s_error$, r_name, 'CANNOT FINISH THE PULSE FILE ' // trim(fel_bank_file_name(run, ihh, 'pulse')))
  return
endif

write (descrip, '(a, i0, a, i0, a)') 'Planes ', nslice + 1, ' to ', nslice + n, &
      ' of the mesh are the banked ones. Each left the window at its own zTransmit and ' // &
      'was drifted to the exit. The planes before them are the window''s.'
call fel_write_bank_group (run%pul_id(ihh), run, ihh, 'pulse', trim(descrip), ferr)
call H5Fclose_f (run%pul_id(ihh), h5e)
run%pul_id(ihh) = 0
if (ferr) return

err_flag = .false.

!------------------------------------------------------------------------------
contains

!+
! Subroutine one_plane (w)
!
! Routine to make w a one-plane wavefront on the window's grid.
!-

subroutine one_plane (w)

type (wavefront_struct) w

!

allocate (w%Ex(size(wfl%Ex, 1), size(wfl%Ex, 2), 1))
if (two) allocate (w%Ey(size(wfl%Ex, 1), size(wfl%Ex, 2), 1))
w%dx = wfl%dx;  w%dy = wfl%dy;  w%dz = wfl%dz
w%wavelength = wfl%wavelength

end subroutine one_plane

!+
! Subroutine read_plane (w, k)
!
! Routine to read banked plane k, in transmission order, into w.
!-

subroutine read_plane (w, k)

type (wavefront_struct) w
integer k

!

call wavefront_openpmd_read_planes (pm_id, 'x', nslice + k, w%Ex, ferr)
if (.not. ferr .and. two) call wavefront_openpmd_read_planes (pm_id, 'y', nslice + k, w%Ey, ferr)

end subroutine read_plane

!+
! Subroutine write_plane (w, i)
!
! Routine to write w as banked plane i in the mesh's order.
!-

subroutine write_plane (w, i)

type (wavefront_struct) w
integer i

!

call wavefront_openpmd_write_planes (pm_id, 'x', nslice + i, w%Ex, ferr)
if (.not. ferr .and. two) call wavefront_openpmd_write_planes (pm_id, 'y', nslice + i, w%Ey, ferr)

end subroutine write_plane

end subroutine fel_finish_pulse

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_write_wake_block (run, z)
!
! Routine to write one block of per-slice wake eloss, z-stamped, to the run's open
! wake file. Written at the hoisted update and at every migration-stride recompute.
! The energy-bookkeeping and stale-wake checks parse these blocks.
!
! Input:
!   run  -- fel_run_struct: Run state. run%coll%wake%eloss holds the current losses.
!   z    -- real(rp): The z position stamped on the block [m].
!
! Output:
!   run  -- fel_run_struct: One block appended to the file on unit run%iu_wake.
!-

subroutine fel_write_wake_block (run, z)

type (fel_run_struct), target :: run

real(rp) z
integer is_w

write (run%iu_wake, '(a, es22.14)') '# z = ', z
do is_w = 1, run%nslice
  write (run%iu_wake, '(i8, es24.16)') is_w, run%coll%wake%eloss(is_w)
enddo

end subroutine fel_write_wake_block

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_write_meta (run, stats_file)
!
! Routine to write the provenance group into the stats file: the resolved input echo
! (every default explicit, straight from the structs), the top-level lattice file's name
! and text, how many files the parser opened, an ISO timestamp, the Lucifer version and
! the Bmad version.
!
! The texts are datasets and not attributes. HDF5 caps a single attribute at 64 kB, and
! the largest that writes is 65495 bytes. A scalar string dataset takes 3 MB. An echoed
! input file runs to 12 kB and a lattice text to 37 kB on a real lattice
! (docs/reading-output.md). A failure would only warn, so the provenance would go missing
! silently on a lattice of ordinary size. Since input_echo carries out_root, two runs
! differ in meta/. The validation suite excludes meta/ by name from identity comparisons.
!
! Note: meta/ is not a reproducibility record. lattice_source is the top-level file
! only. Bmad's "call, file =" pulls in more, and a wrapper lattice records a call
! statement without the lattice it calls. n_lattice_files says how many files the parser
! opened, so a reader can see whether the text is the whole story. For reproduction the
! file offers the lattice/ table (every tracked element with the values the physics
! used) beside the input echo.
!
! Serializing the lattice is not an option. write_bmad_lattice_file inlines a grid_field
! as ASCII under one_file$ and writes sibling binary files otherwise. So no output_form
! is both complete and bounded, and Bmad has no HDF5 lattice format. lat%creation_hash is
! not a substitute either. It hashes inode and size, which makes it machine specific and
! not a fingerprint of content.
!
! A stats file is meant to travel, so nothing here identifies a person by default. The
! timestamp and the two versions identify the run. The login name and working directory
! go in only under run.record_environment. The file records a lattice basename.
! The input file is echoed as read, so an absolute path there is still an absolute path
! in the file. Genesis always records user and cwd. Parity is not a reason to leak.
!
! Note: This is best effort. A failure warns, and never fails the run.
!
! Input:
!   run         -- fel_run_struct: Run state.
!   stats_file  -- character(*): The stats file to reopen and annotate.
!
! Output:
!   None. Failures are reported with a warning and otherwise ignored.
!-

subroutine fel_write_meta (run, stats_file)

type (fel_run_struct), target :: run
character(*) stats_file
character(:), allocatable :: txt
character(24) stamp
character(200) user_name, cwd
character(400) base, path
character(8) date_s
character(10) time_s
integer(hid_t) f_id, g_id
integer h5e, nfile, ix
logical merr
character(*), parameter :: r_name = 'fel_write_meta'

!

call hdf5_open_file (stats_file, 'APPEND', f_id, merr)
if (merr) then
  call out_io (s_warn$, r_name, 'Could not reopen the stats file for meta/.')
  return
endif
call H5Gcreate_f (f_id, 'meta', g_id, h5e)
if (h5e < 0) then
  call H5Fclose_f (f_id, h5e)
  return
endif
call H5LTset_attribute_string_f (f_id, 'meta', 'kind', 'provenance', h5e)
call H5LTset_attribute_string_f (f_id, 'meta', 'description', &
      'How this file was made. Which lattice, not the lattice itself.', h5e)

merr = .false.
call resolved_input_text (run, txt)
call fel_h5_text (g_id, 'input_echo', 'input echo', &
      'The run''s input with every default made explicit, written from the structs. ' // &
      'File names appear as the user typed them.', txt, merr)

! The lattice: which one, and how much of it is here. bmad_parser opens the top-level
! file plus every file its "call, file =" statements reach. file_text reads only the
! top-level file, so a wrapper lattice records a call statement and drops what it calls. The
! count says whether the text is complete. It comes from the parser's own tally, which
! outlives the parse because parser_end_stuff deallocates only the array of names. Zero
! means the tally was not available, and not that no file was read.

ix = splitfilename(run%lat_file, path, base)
call fel_h5_text (g_id, 'lattice_file', 'lattice file', &
      'Base name of the top-level lattice file. The directory is deliberately absent: ' // &
      'a stats file travels, and a path is machine local.', base, merr)

call file_text (trim(run%lat_file), txt)
call fel_h5_text (g_id, 'lattice_source', 'lattice source', &
      'Text of the top-level lattice file only. A file reached by "call, file =" is ' // &
      'NOT here, so this is not a reproducibility record. See n_lattice_files, and ' // &
      'read lattice/ for the tracked line as the physics used it.', txt, merr)

nfile = max(bp_com%num_lat_files, 0)
call fel_h5_int (g_id, 'n_lattice_files', '1', 'lattice files', &
      'How many files the parser opened, so one means lattice_source is the whole ' // &
      'story. Zero means the parser did not report a count.', '', nfile, merr)

call date_and_time (date_s, time_s)
stamp = date_s(1:4) // '-' // date_s(5:6) // '-' // date_s(7:8) // 'T' // &
        time_s(1:2) // ':' // time_s(3:4) // ':' // time_s(5:6)
call fel_h5_text (g_id, 'timestamp', 'timestamp', &
      'Local time the run finished, ISO 8601 without a zone.', stamp, merr)
call fel_h5_text (g_id, 'lucifer_version', 'Lucifer version', &
      'The Lucifer version the run was built from: the release tag, git describe of the ' // &
      'commit, or what the build was told. A version ending in -dirty was built from a ' // &
      'tree with uncommitted changes.', lucifer_version$, merr)
call fel_h5_int (g_id, 'bmad_inc_version', '1', 'Bmad lattice-format version', &
      'The bmad_inc_version$ the run was built against.', '', bmad_inc_version$, merr)
call fel_h5_text (g_id, 'bmad_version', 'Bmad distribution version', &
      'The Bmad distribution version the run was built from, the string Tao reports. ' // &
      'The date-like form is Bmad''s own, and bmad_inc_version above is a different ' // &
      'number: the lattice-format version.', bmad_version_date, merr)

! The environment, only when asked for. These identify a person and a machine, and they
! are of no use to a reader elsewhere.

if (run%global%record_environment) then
  call get_environment_variable ('USER', user_name)
  call fel_h5_text (g_id, 'user', 'user', &
        'Who ran it, from $USER. Present because run.record_environment was set.', &
        user_name, merr)
  call get_environment_variable ('PWD', cwd)
  call fel_h5_text (g_id, 'cwd', 'directory', &
        'Where it ran, from $PWD. Present because run.record_environment was set.', &
        cwd, merr)
endif

if (merr) call out_io (s_warn$, r_name, 'Could not write all of meta/.')

call H5Gclose_f (g_id, h5e)
call H5Fclose_f (f_id, h5e)

end subroutine fel_write_meta

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_write_params (run, stats_file)
!
! Routine to write the input tree into the stats file: params/, one subgroup per input
! structure the program honors. Each carries @struct and holds the resolved value of
! every honored component. So the inputs are data, two runs diff by their inputs, and
! no reader needs a defaults table. What the run produced is run/, written by
! fel_stats_write. The literal text stays in meta/input_echo.
!
! Components the program refuses are not here, since recording a knob that did nothing
! would misstate the run. The quiet start's honored beam_init set is the table in
! docs/input-reference.md. List-valued inputs (the dump locators, the harmonic set, the
! field files) are string or integer array attributes on their subgroup. A dataset
! needs an axis, and a list of inputs has none worth defining. A blank file-name input
! is omitted: absence means unset. bmad_com and space_charge_com are written whole,
! since the namelist exposes them whole and all of Bmad's tracking reads them.
!
! Note: This is best effort, like fel_write_meta. A failure warns, and never fails the run.
!
! Input:
!   run         -- fel_run_struct: Run state, after setup.
!   stats_file  -- character(*): The stats file to reopen and annotate.
!-

subroutine fel_write_params (run, stats_file)

type (fel_run_struct), target :: run
type (fel_key_struct), allocatable :: keys(:)
character(*) stats_file
integer(hid_t) f_id, p_id, g_id, b_id
integer h5e, i, it
logical merr, there
character(16) table
character(12), parameter :: origin_name(0:3) = [character(12):: 'default', 'input file', 'command line', 'derived']
character(*), parameter :: r_name = 'fel_write_params'

!

call hdf5_open_file (stats_file, 'APPEND', f_id, merr)
if (merr) then
  call out_io (s_warn$, r_name, 'Could not reopen the stats file for params/.')
  return
endif
call H5Gcreate_f (f_id, 'params', p_id, h5e)
if (h5e < 0) then
  call H5Fclose_f (f_id, h5e)
  return
endif
merr = .false.
call H5LTset_attribute_string_f (f_id, 'params', 'kind', 'input', h5e)
call H5LTset_attribute_string_f (f_id, 'params', 'description', &
      'What the user controlled: one subgroup per table of the input file, every key ' // &
      'resolved after defaults, from the schema (schema version ' // fel_schema_version$ // ').', h5e)

call fel_schema_table (keys)

! One subgroup per table, [beam.resample] inside beam/. Every recorded key is written by
! its kind: a number or a flag as a scalar dataset, a word as a string dataset, a path as
! text when it is set, and a list as an attribute of the subgroup, since a list of inputs
! has no axis worth defining. out_root is the one key not recorded: it is the run's own
! name and not a physics input, every output file carries it, and as data it would make
! two otherwise identical runs compare different.

do it = 1, size(fel_schema_tables$)
  table = fel_schema_tables$(it)
  ! Each subgroup carries @kind and @description, as the stats specification requires of
  ! a group.
  if (table == 'beam.resample') then
    call H5Gopen_f (p_id, 'beam', b_id, h5e)
    call H5Gcreate_f (b_id, 'resample', g_id, h5e)
    call H5LTset_attribute_string_f (b_id, 'resample', 'kind', 'input', h5e)
    call H5LTset_attribute_string_f (b_id, 'resample', 'description', 'The input file''s [beam.resample] table.', h5e)
  else
    call H5Gcreate_f (p_id, trim(table), g_id, h5e)
    call H5LTset_attribute_string_f (p_id, trim(table), 'kind', 'input', h5e)
    call H5LTset_attribute_string_f (p_id, trim(table), 'description', &
          'The input file''s [' // trim(table) // '] table.', h5e)
  endif
  merr = merr .or. (h5e < 0)
  do i = 1, size(keys)
    if (keys(i)%table /= table .or. .not. keys(i)%recorded) cycle
    call write_key (keys(i))
    ! Where the value came from, on every key that is a dataset.
    call H5Lexists_f (g_id, trim(keys(i)%name), there, h5e)
    if (there .and. allocated(run%key_origin)) call H5LTset_attribute_string_f (g_id, trim(keys(i)%name), &
          'origin', trim(origin_name(run%key_origin(i))), h5e)
  enddo
  call H5Gclose_f (g_id, h5e)
  if (table == 'beam.resample') call H5Gclose_f (b_id, h5e)
enddo

if (merr) call out_io (s_warn$, r_name, 'Could not write all of params/.')
call H5Gclose_f (p_id, h5e)
call H5Fclose_f (f_id, h5e)

!-------------------------------------------------------
contains

subroutine write_key (key)

type (fel_key_struct) key
real(rp) rv
integer iv
logical lv
character(16) unit
character(:), allocatable :: sv
integer, allocatable :: ilist(:)
character(400), allocatable :: slist(:)

!

call fel_schema_value (run, key, rv, iv, lv, sv, ilist, slist)

select case (key%kind)
case (key_real$)
  ! A dimensionless key carries the unit 1, which the stats specification requires.
  unit = key%unit
  if (unit == '') unit = '1'
  call fel_h5_real (g_id, trim(key%name), trim(unit), trim(key%label), trim(key%descrip), '', rv, merr)
case (key_int$)
  call fel_h5_int (g_id, trim(key%name), '1', trim(key%label), trim(key%descrip), '', iv, merr)
case (key_logical$)
  call fel_h5_flag (g_id, trim(key%name), trim(key%label), trim(key%descrip), lv, merr)
case (key_string$)
  call fel_h5_str (g_id, trim(key%name), trim(key%label), trim(key%descrip), '', [sv], merr)
case (key_path$)
  if (sv /= '') call fel_h5_text (g_id, trim(key%name), trim(key%label), trim(key%descrip) // &
        ' As the input file wrote it.', sv, merr)
case (key_int_list$)
  if (size(ilist) > 0) then
    call H5LTset_attribute_int_f (g_id, '.', trim(key%name), ilist, int(size(ilist), size_t), h5e)
    merr = merr .or. (h5e < 0)
  endif
case (key_string_list$)
  ! One text, an entry a line, since a path can be longer than a string dataset holds.
  if (size(slist) > 0) then
    sv = trim(slist(1))
    do iv = 2, size(slist)
      sv = sv // char(10) // trim(slist(iv))
    enddo
    call fel_h5_text (g_id, trim(key%name), trim(key%label), trim(key%descrip) // ' One entry a line.', sv, merr)
  endif
end select

end subroutine write_key

end subroutine fel_write_params


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_write_lattice (run, stats_file)
!
! Routine to write the lattice table into the stats file: one row per tracked element,
! which is what a layout plot needs. Genesis writes its Lattice/ group as per-step
! arrays. This is a table instead, joined to the records through coords/ix_ele. So
! nothing is duplicated per record, and the two cannot drift.
!
! Note: This is *not* a lattice serialization. It carries what a plot and a join need.
! Anything that needs the real lattice reads meta/lattice_file or meta/lattice_source.
! Bmad has no portable HDF5 lattice format, and a statistics file is the wrong place to
! invent one.
!
! The table's axis is named ele and gets its own coordinate, coords/ele = 0..n_ele.
! coords/ix_ele is then a variable on the record axis whose values index that one. That
! is the join. One name for two axes of different length would be a collision every
! reader resolves by hand.
!
! aw is the value the physics used, derived through b_max and l_period with the helical
! or planar factor, and not the raw attribute. That is the number a Genesis user
! compares against Lattice/aw. It is taken from the run, since deriving it here would
! repeat fel_setup_lattice.
!
! Note: This is best effort, like fel_write_meta. A failure warns, and never fails the run.
!
! Input:
!   run         -- fel_run_struct: Run state, after setup.
!   stats_file  -- character(*): The stats file to reopen and annotate.
!-

subroutine fel_write_lattice (run, stats_file)

type (fel_run_struct), target :: run
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
character(*) stats_file
character(fel_h5_str_len$), allocatable :: names(:), keys(:), tmethod(:), scmethod(:)   ! 0:n_ele_track
real(rp), allocatable :: s1(:), s2(:), len_(:), dstep(:), b_max(:), aw(:), l_per(:)
real(rp), allocatable :: ku(:), k1(:), tilt(:), z_off(:)
integer, allocatable :: is_fel(:), helical(:)
integer(hid_t) f_id, g_id
integer h5e, ie, ne
logical merr
character(*), parameter :: r_name = 'fel_write_lattice'

!

! The table is indexed by ix_ele, element 0 (Bmad's beginning element) included. So
! coords/ix_ele indexes it with no offset arithmetic: lattice/name(ix_ele) is the
! record's element. The first record of a run sits at element 0.

branch => run%lat%branch(0)
ne = branch%n_ele_track

allocate (names(0:ne), keys(0:ne), tmethod(0:ne), scmethod(0:ne), s1(0:ne), s2(0:ne), len_(0:ne), dstep(0:ne))
allocate (b_max(0:ne), aw(0:ne), l_per(0:ne), ku(0:ne), k1(0:ne), tilt(0:ne), z_off(0:ne))
allocate (is_fel(0:ne), helical(0:ne))

do ie = 0, ne
  ele => branch%ele(ie)
  names(ie) = ele%name
  keys(ie) = key_name(ele%key)
  tmethod(ie) = tracking_method_name(ele%tracking_method)
  scmethod(ie) = space_charge_method_name(ele%space_charge_method)
  s1(ie) = ele%s_start
  s2(ie) = ele%s
  len_(ie) = ele%value(l$)
  dstep(ie) = ele%value(ds_step$)
  b_max(ie) = ele%value(b_max$)
  l_per(ie) = ele%value(l_period$)
  k1(ie) = ele%value(k1$)
  tilt(ie) = ele%value(tilt_tot$)
  z_off(ie) = ele%value(z_offset_tot$)
  is_fel(ie) = 0
  if (ie > 0) is_fel(ie) = merge(1, 0, run%is_fel(ie))
  aw(ie) = 0;  ku(ie) = 0;  helical(ie) = 0
  if (is_fel(ie) == 1) then
    aw(ie) = run%und_of(ie)%aw
    ku(ie) = run%und_of(ie)%ku
    helical(ie) = merge(1, 0, run%und_of(ie)%helical)
  endif
enddo

call hdf5_open_file (stats_file, 'APPEND', f_id, merr)
if (merr) then
  call out_io (s_warn$, r_name, 'Could not reopen the stats file for lattice/.')
  return
endif
call H5Gcreate_f (f_id, 'lattice', g_id, h5e)
if (h5e < 0) then
  call H5Fclose_f (f_id, h5e)
  return
endif

merr = .false.
call H5LTset_attribute_string_f (f_id, 'lattice', 'kind', 'table', h5e)
call H5LTset_attribute_string_f (f_id, 'lattice', 'description', &
      'One row per tracked element, indexed by the ele axis.', h5e)

call fel_h5_str (g_id, 'name', 'name', &
      'Element name. Indexed BY the ele axis, element 0 included.', 'ele', names, merr)
call fel_h5_str (g_id, 'key', 'class', 'Bmad element class, as key_name gives it.', &
      'ele', keys, merr)
call fel_h5_real (g_id, 's_start', 'm', 's start', 'Upstream end of the element.', &
      'ele', s1, merr)
call fel_h5_real (g_id, 's_end', 'm', 's end', 'Downstream end of the element.', 'ele', s2, merr)
call fel_h5_real (g_id, 'l', 'm', 'length', 'Element length.', 'ele', len_, merr)
call fel_h5_real (g_id, 'ds_step', 'm', 'ds_step', &
      'Integration step the walk used, which is what sets the record density.', &
      'ele', dstep, merr)
call fel_h5_flag (g_id, 'is_fel', 'is FEL', &
      'One where the element is an FEL segment the FEL step tracked.', 'ele', is_fel, merr)
call fel_h5_str (g_id, 'tracking_method', 'tracking method', &
      'How the element was tracked, as Bmad names it. FEL_Averaged and FEL_Unaveraged ' // &
      'are the FEL segments, and every other element carries the method Bmad used.', &
      'ele', tmethod, merr)
call fel_h5_str (g_id, 'space_charge_method', 'space charge method', &
      'The element''s own space_charge_method, as Bmad names it, which is what decided ' // &
      'whether space charge acted here. Off unless the element asked for it, and ' // &
      'Bmad''s csr_and_space_charge_on had to be true as well.', 'ele', scmethod, merr)
call fel_h5_real (g_id, 'b_max', 'T', 'b_max', 'Peak undulator field, zero elsewhere.', &
      'ele', b_max, merr)
call fel_h5_real (g_id, 'aw', '1', 'aw', &
      'Rms undulator parameter as the physics used it: c*b_max/(ku*m_e c^2), ' // &
      'divided by sqrt(2) for a planar device. Zero off an FEL segment.', 'ele', aw, merr)
call fel_h5_real (g_id, 'l_period', 'm', 'period', 'Undulator period.', 'ele', l_per, merr)
call fel_h5_real (g_id, 'ku', '1/m', 'ku', 'Undulator wavenumber, twopi/l_period.', &
      'ele', ku, merr)
call fel_h5_flag (g_id, 'helical', 'helical', &
      'One for a helical device, zero for planar.', 'ele', helical, merr)
call fel_h5_real (g_id, 'k1', '1/m^2', 'k1', &
      'Quadrupole strength, signed as Bmad signs it.', 'ele', k1, merr)
call fel_h5_real (g_id, 'tilt', 'rad', 'tilt', &
      'Element tilt. On a planar FEL segment this is the wiggle-plane rotation, ' // &
      'which is the polarization spec.', 'ele', tilt, merr)
call fel_h5_real (g_id, 'z_offset', 'm', 'z offset', &
      'Longitudinal misalignment, the inter-segment phasing knob.', 'ele', z_off, merr)
call H5Gclose_f (g_id, h5e)

! The axis itself, into the coords/ group. Every name in an @axes attribute resolves to
! a coords/ dataset, and only this routine knows the length of the lattice/ axis.

call H5Gopen_f (f_id, 'coords', g_id, h5e)
if (h5e < 0) then
  merr = .true.
else
  call fel_h5_int (g_id, 'ele', '1', 'element', &
        'Lattice element index, the axis the lattice/ table runs over. coords/ix_ele ' // &
        'is a variable whose values index it.', 'ele', [(ie, ie = 0, ne)], merr)
  call H5Gclose_f (g_id, h5e)
endif

if (merr) call out_io (s_warn$, r_name, 'Could not write all of lattice/.')

call H5Fclose_f (f_id, h5e)

end subroutine fel_write_lattice

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine resolved_input_text (run, txt)
!
! Routine to render the resolved input file as one string, every key of the schema with its
! value, through fel_schema_write_toml and a scratch file.
!
! Input:
!   run  -- fel_run_struct: Run state.
!
! Output:
!   txt  -- character(:), allocatable: The resolved input echo.
!-

subroutine resolved_input_text (run, txt)

type (fel_run_struct), target :: run
character(:), allocatable :: txt
integer iu

open (newunit = iu, status = 'scratch', action = 'readwrite')
call fel_schema_write_toml (run, iu)
call unit_text (iu, txt)
close (iu)

end subroutine resolved_input_text

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine file_text (file_name, txt)
!
! Routine to read a whole file into one string. An unreadable file yields the
! string '(unreadable)'.
!
! Input:
!   file_name  -- character(*): File to read.
!
! Output:
!   txt        -- character(:), allocatable: The file's text.
!-

subroutine file_text (file_name, txt)

character(*) file_name
character(:), allocatable :: txt
integer iu, ios

open (newunit = iu, file = file_name, status = 'old', action = 'read', iostat = ios)
if (ios /= 0) then
  txt = '(unreadable)'
  return
endif
call unit_text (iu, txt)
close (iu)

end subroutine file_text

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine unit_text (iu, txt)
!
! Routine to read an open unit from the start into one string.
!
! Input:
!   iu   -- integer: Open unit. Rewound before reading.
!
! Output:
!   txt  -- character(:), allocatable: The unit's text.
!-

subroutine unit_text (iu, txt)

integer iu, ios
character(:), allocatable :: txt
character(1000) line

rewind (iu)
txt = ''
do
  read (iu, '(a)', iostat = ios) line
  if (ios /= 0) exit
  txt = txt // trim(line) // new_line('a')
enddo

end subroutine unit_text

end module fel_io_mod
