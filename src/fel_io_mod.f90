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
use wavefront_openpmd_mod
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
! For bp_com%num_lat_files, the parser's tally of how many files it opened. Read only in
! fel_write_meta.
use bmad_parser_struct, only: bp_com
! For bmad_version_date, the distribution's own version. Tao reports the same string.
! Lucifer builds against the same distribution, so there is no second version to keep.
use bmad_version_mod, only: bmad_version_date
!$ use omp_lib

implicit none

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
! The field's records are rotated to time order to be written and rotated back, so the
! run continues from the state it had. fel_dump_field_set leaves them unrotated. That is
! right at the end of a run and wrong in the middle of one. A diagnostic that changes
! the state it observes is the failure the FP32 twin's read-only proof exists to catch.
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
type (wavefront_struct) wf_sub
complex(wf_rp), allocatable :: eslab(:,:,:)
real(rp) s_frame
integer ie, ihh, first_was, ifr
logical at_end, err_flag, eerr, whole, ckpt
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

    ! A scratch copy. Inside an unaveraged segment the record is the substep's midpoint
    ! field, half a substep past the plane this frame is named for (fel_unaveraged_mod).
    ! The frame carries the plane's own field. So the copy is carried back and written,
    ! and the live record is not touched. A run with frames on and one with them off then
    ! hold the same field to the bit. Carrying the live record back and forward again
    ! would leave two transform pairs of rounding in it, 5e-15 of the field on the mixed
    ! line (docs/validation.md). A slice range takes its slab from the same copy. The copy
    ! is in time order for the write, as the whole-window path below rotates the live
    ! record to be.

    wf_sub = ffield(ihh)%wf
    if (run%ustate%active) then
      do ifr = 1, size(wf_sub%Ex, 3)
        call fel_field_diffract_by (wf_sub, ifr, -0.5_rp * run%ustate%dsub, eerr)
        if (eerr) return
      enddo
    endif
    if (first_was /= 0) then
      wf_sub%Ex = cshift(wf_sub%Ex, shift = first_was, dim = 3)
      if (allocated(wf_sub%Ey)) wf_sub%Ey = cshift(wf_sub%Ey, shift = first_was, dim = 3)
    endif
    if (.not. whole) then

      ! The slab keeps its place. The mesh's z origin is the first slice's. So a reader
      ! that walks the standard attributes puts the cut field where the whole window's
      ! would sit. It sits beside the cut beam, whose timeOffset places each slice.

      eslab = wf_sub%Ex(:, :, run%dump_is1:run%dump_is2)
      call move_alloc (eslab, wf_sub%Ex)
      if (allocated(wf_sub%Ey)) then
        eslab = wf_sub%Ey(:, :, run%dump_is1:run%dump_is2)
        call move_alloc (eslab, wf_sub%Ey)
      endif
    endif
    call wavefront_write_openpmd (wf_sub, trim(prefix) // trim(hsuf) // '.wf.h5', eerr, &
                                  z_offset = (run%dump_is1 - 1) * ffield(ihh)%wf%dz)
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
character(200) line
character(*), parameter :: r_name = 'lucifer'
character(*), parameter :: rule = &
      '================================================================================'
integer nfel, n_omp

!

branch => run%lat%branch(0)
nfel = count(run%is_fel)
n_omp = 1
!$ n_omp = omp_get_max_threads()

call out_io (s_blank$, r_name, rule)
call out_io (s_blank$, r_name, ' Lucifer -- FEL tracking in Bmad, Bmad version ' // bmad_version_date)
call out_io (s_blank$, r_name, &
      '--------------------------------------------------------------------------------')

write (line, '(a, a)') ' Lattice     ', trim(run%lat_file)
call out_io (s_blank$, r_name, trim(line))
write (line, '(a, i0, a, f0.3, a, i0, a)') '             ', branch%n_ele_track, ' elements, ', &
      branch%ele(branch%n_ele_track)%s, ' m, ', nfel, ' FEL segments'
call out_io (s_blank$, r_name, trim(line))

if (run%nslice == 1) then
  write (line, '(a, i0, a, f0.2)') ' Beam        1 slice x ', run%fbeam%slice(1)%n, &
        ' particles, gamma0 = ', run%gamma0
elseif (all(run%fbeam%slice%n == run%fbeam%slice(1)%n)) then
  write (line, '(a, i0, a, i0, a, f0.2)') ' Beam        ', run%nslice, ' slices x ', &
        run%fbeam%slice(1)%n, ' particles per slice, gamma0 = ', run%gamma0
else
  write (line, '(a, i0, a, i0, a, i0, a, i0, a, f0.2)') ' Beam        ', run%nslice, ' slices, ', &
        sum(run%fbeam%slice%n), ' particles, ', minval(run%fbeam%slice%n), ' to ', &
        maxval(run%fbeam%slice%n), ' per slice, gamma0 = ', run%gamma0
endif
call out_io (s_blank$, r_name, trim(line))

! The grid is the built wavefront's and not the input's. A field read from a file
! carries its own grid, and the deck then states none.

block
  integer n_grid(3)
  n_grid = wavefront_shape(run%ffield(1)%wf)
  write (line, '(5a, i0, a, i0, 2a)') ' Radiation   lambda0 = ', &
        trim(adjustl(fel_si_str(run%winit%lambda0, 'm'))), &
        ', slice spacing ', trim(adjustl(fel_si_str(run%fbeam%slice_spacing, 'm'))), ', ', &
        run%n_harm, ' field(s), grid ', n_grid(1), ' points of ', &
        trim(adjustl(fel_si_str(run%ffield(1)%wf%dx, 'm')))
end block
call out_io (s_blank$, r_name, trim(line))

! What the run will do and how long it will take. The count is exact. The estimate is
! one machine's, which the line says (fel_cost_mod). load_only stops after this block,
! so a deck's cost can be seen without tracking it (docs/user-guide.md).

block
  real(rp) n_pstep, n_gpt, secs
  integer nstep, n_particle, ngrid, nstep_unavg, nbeamlet, nharm_noise
  character(16) backend
  character(64) ramp_str

  call fel_cost_count (run, nstep, n_particle, n_pstep, n_gpt, ngrid, nstep_unavg)
  write (line, '(a, i0, a, i0, a, es8.2, a, es8.2, a)') ' Work        ', nstep, &
        ' FEL steps, ', n_particle, ' macroparticles, ', n_pstep, ' particle-steps, ', &
        n_gpt, ' grid points per step'
  call out_io (s_blank$, r_name, trim(line))

  ! The device's transform pair is single precision. It loses a fixed fraction of the
  ! field's energy every time it runs, always in one direction. So what a deck loses is
  ! set by how often the pair runs: once a record step in an averaged element, and nsub
  ! times in an unaveraged one (nstep already counts the substeps). Nothing else in the
  ! output says which regime a deck is in (validation.md, val-device-unaveraged).

  if (run%dev%on) then
    write (line, '(a, i0, a, i0, a, es8.2, a, es8.2, a)') ' Device      FP32 transform pair ', nstep, &
          ' times a field record (', nstep_unavg, ' unaveraged substeps), projected field-energy loss ', &
          fel_dev_pair_loss(ngrid), ' a pair, ', fel_dev_pair_loss(ngrid) * nstep, ' over the line'
    call out_io (s_blank$, r_name, trim(line))
  endif

  ! The interlude's pieces, only when run.interlude_ds_step cuts them. A run that
  ! tracks an interlude whole prints nothing here.

  if (run%n_int_ele > 0) then
    write (line, '(3a, i0, a, i0, a)') ' Interlude   steps of ', &
          trim(adjustl(fel_si_str(run%global%interlude_ds_step, 'm'))), ': ', &
          run%n_int_piece, ' pieces, ', run%n_int_ele, ' element(s) cut'
    call out_io (s_blank$, r_name, trim(line))
  endif

  ! The unaveraged elements' substeps and ramps, only when the line has one. The ramp
  ! says where it came from, since a hard edge brings back the handoff's phase jump.

  if (any(run%fel_mode == unaveraged$)) then
    select case (run%global%unaveraged_ramp_periods)
    case (-1)
      ramp_str = 'a hard edge, the sentinel -1'
    case (fel_ramp_periods_default$)
      write (ramp_str, '(a, i0, a)') 'sin^2 ramps of ', run%global%unaveraged_ramp_periods, &
            ' periods at each end, the default'
    case default
      write (ramp_str, '(a, i0, a)') 'sin^2 ramps of ', run%global%unaveraged_ramp_periods, &
            ' periods at each end, as the deck sets'
    end select
    write (line, '(a, i0, 2a)') ' Unaveraged  ', run%global%unaveraged_steps_per_period, &
          ' substeps a period, ', trim(ramp_str)
    call out_io (s_blank$, r_name, trim(line))
  endif

  call fel_cost_estimate (run, secs, backend)
  write (line, '(5a)') ' Estimate    ', trim(fel_cost_secs_str(secs)), ' of walk on the ' // &
        trim(backend) // ' path, for ', fel_cost_machine$, '.'
  call out_io (s_blank$, r_name, trim(line))

  ! The frame series, when it is on. A comb of zero on a long line is a frame per step.
  ! So the bytes are stated before the run rather than discovered on a full disk. The
  ! per-frame size is the beam's coordinates and the field set's grid, the two counts
  ! above. Only the beam and the field carry the slice range, so their sizes scale with
  ! it. The reductions are per record whether or not frames are on, and the transverse
  ! projection is the grid's own size. So they are sized separately, and stated even
  ! without a series.

  block
    real(rp) n_rng, b_frame, f_frame, red_bytes
    integer nrng
    character(12) lab
    nrng = run%dump_is2 - run%dump_is1 + 1
    n_rng = real(nrng, rp) / max(1, run%nslice)
    ! 68 bytes a macroparticle (docs/performance.md): six coordinates, the time the phase
    ! folds into, the offset that places the slice, the reference momentum and the
    ! label. The writer stores a uniform weight and the zero z as constant records, which
    ! cost nothing per particle.
    b_frame = 68.0_rp * n_particle * n_rng
    f_frame = 16.0_rp * n_gpt * n_rng
    red_bytes = 0
    if (run%stats%n_red > 0) red_bytes = 8.0_rp * run%stats%nrec * run%stats%n_red * &
          (real(ngrid, rp)**2 + 2.0_rp * run%nslice * ngrid)

    if (run%global%dump_at_comb) then
      write (line, '(a, i0, a, i0, a, es8.2, a)') ' Frames      ', run%stats%nrec, &
            ' at the comb over ', nrng, ' slices, about ', &
            run%stats%nrec * (b_frame + f_frame), ' bytes of beam and field'
      call out_io (s_blank$, r_name, trim(line))
    endif
    if (red_bytes > 0) then
      lab = ' Reduced    '
      if (run%global%dump_at_comb) lab = '            '
      write (line, '(a, i0, a, es8.2, a)') lab // ' The field''s reductions over ', &
            run%stats%n_red, ' of the field set, about ', red_bytes, &
            ' bytes in the stats file'
      call out_io (s_blank$, r_name, trim(line))
    endif
  end block

  ! The load against the floor of the filter's convergence: 1024 macroparticles per
  ! slice in 128 beamlets. Below it, the mode power in the exponential regime rises on
  ! its own (docs/fel-physics.md sec-convergence).

  nbeamlet = minval(run%fbeam%slice%n) / max(1, run%bparam%beamlet_size)
  write (line, '(a, i0, a, i0, a, i0, a)') ' Load        ', minval(run%fbeam%slice%n), &
        ' macroparticles in the thinnest slice, beamlets of ', run%bparam%beamlet_size, &
        ', ', nbeamlet, ' beamlets.'
  call out_io (s_blank$, r_name, trim(line))

  if (.not. run%global%source_filter) then
    call out_io (s_blank$, r_name, '             The source filter is off, and no estimate of the ' // &
          'wide-angle share carries')
    call out_io (s_blank$, r_name, '             between machines. The two tests are a run at twice ' // &
          'this load and the')
    call out_io (s_blank$, r_name, '             convergence report in this run''s footer.')
  else if (minval(run%fbeam%slice%n) >= 1024 .and. nbeamlet >= 128) then
    call out_io (s_blank$, r_name, '             This load is the one the filter''s convergence ' // &
          'was measured on (fel-physics.md).')
  else
    call out_io (s_blank$, r_name, '             Under the measured load of 1024 macroparticles in ' // &
          '128 beamlets, where the')
    call out_io (s_blank$, r_name, '             exponential-regime mode power rises on its own and ' // &
          'the filter is not')
    call out_io (s_blank$, r_name, '             measured to remove it (fel-physics.md).')
  endif

  ! The beamlet size, stated by what it does. The quiet start leaves no bunching below
  ! harmonic beamlet_size. fel_fawley_noise then imposes the physical level on harmonics
  ! 1 through (beamlet_size - 1)/2.

  nharm_noise = (run%bparam%beamlet_size - 1) / 2
  if (run%bparam%shot_noise) then
    write (line, '(a, i0, a)') '             Shot noise on harmonics 1 to ', nharm_noise, &
          ', which is (beamlet_size - 1)/2.'
  else
    write (line, '(a, i0, a)') '             No shot noise. The quiet start leaves no bunching ' // &
          'below harmonic ', run%bparam%beamlet_size, '.'
  endif
  call out_io (s_blank$, r_name, trim(line))
end block

write (line, '(a, l1, a, l1, a, l1)') ' Switches    sr wakes ', run%coll%wake%on, &
      ', space charge ', any(run%sc_here), ', radiation damping ', bmad_com%radiation_damping_on
call out_io (s_blank$, r_name, trim(line))

write (line, '(3a, i0)') ' Output      out_root "', trim(run%global%out_root), '", threads ', n_omp
call out_io (s_blank$, r_name, trim(line))
call out_io (s_blank$, r_name, rule)

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
! Input:
!   run -- fel_run_struct: The finished run state.
!-

subroutine fel_write_footer (run)

type (fel_run_struct), target :: run
type (fel_stats_struct), pointer :: stats
character(400) line
character(*), parameter :: r_name = 'lucifer'
type (fel_convergence_struct) cvg
integer ir, ih, n_listed, is, n_thin, n_shown, nbl
real(rp) pow, ene, bun, worst, secs, ratio_s
character(8) hsuf
character(16) backend
integer, parameter :: n_thin_max$ = 8

!

stats => run%stats
call out_io (s_blank$, r_name, &
      '--------------------------------------------------------------------------------')

write (line, '(a, f0.3, a, i0, a)') ' Done        ', run%z_now, ' m, ', run%stats%iend, ' element ends'
call out_io (s_blank$, r_name, trim(line))

! The exit light, from whichever row holds the run's last state.

call fel_stats_exit_light (stats, pow, ene, bun)
write (line, '(5a, f6.4)') ' Exit        power ', trim(adjustl(fel_si_str(pow, 'W'))), &
      ', pulse energy ', trim(adjustl(fel_si_str(ene, 'J'))), ', <|b|> ', bun
call out_io (s_blank$, r_name, trim(line))

! The mode against the wide-angle emission of the point beamlets, from the records the
! run already took (docs/startup-noise.md). The share is reported and not estimated. A
! formula fitted on one machine is 40 times low on another (docs/fel-physics.md), while
! the split itself holds on all three.

call fel_stats_convergence (stats, cvg)
if (cvg%ok) then
  write (line, '(4a)') ' Split       ', trim(adjustl(fel_si_str(cvg%angle, 'rad'))), &
        ', ', trim(run%split_origin) // '. Inside it the mode, outside it the beamlets.'
  call out_io (s_blank$, r_name, trim(line))
  do ir = 1, 2
    if (ir == 2 .and. cvg%z(2) == cvg%z(1)) cycle
    write (line, '(2a, f0.3, 5a, es9.2, a, f6.4)') '             ', &
          merge('mode peak, z =   ', 'last record, z = ', ir == 1), cvg%z(ir), ' m: inside ', &
          trim(adjustl(fel_si_str(cvg%p_in(ir), 'W'))), ', outside ', &
          trim(adjustl(fel_si_str(cvg%p_out(ir), 'W'))), ', ratio ', cvg%ratio(ir), &
          ', <|b|> ', cvg%bunching(ir)
    call out_io (s_blank$, r_name, trim(line))
  enddo

  ! The verdict reads the worse of the two records. The criterion is stated at the mode's
  ! saturation. A run that passes there can still exit with most of its power outside
  ! the mode, since the artifact keeps growing after the mode turns over.

  worst = maxval(cvg%ratio)
  if (worst > 0.1_rp) then
    call out_io (s_warn$, r_name, &
          'The power outside the split angle reaches \f0.2\ times the power inside it.', &
          'That part of the total is the wide-angle emission of the point beamlets, which', &
          'no physical beam radiates, so a total power quoted from this run is mostly', &
          'representation. Possible solution: more macroparticles per slice at the same', &
          'beamlet_size. The criterion and the measured levels are in docs/startup-noise.md.', &
          r_array = [worst])
    if (.not. run%global%source_filter) then
      call out_io (s_warn$, r_name, 'This run has run.source_filter = false. The filter is on by ' // &
            'default and removes that emission at a converged load.')
    endif
  else if (worst < 0.01_rp) then
    call out_io (s_info$, r_name, 'The power outside the split angle stays under \f0.4\ ' // &
          'times the power inside it.', &
          'The total power is converged by the criterion of docs/startup-noise.md.', &
          r_array = [worst])
  endif

  ! The thin slices, named rather than averaged away. A window ratio is a mean over the
  ! window. A tail under the floor of 128 beamlets, or whose light is mostly wide angle,
  ! disappears into it. This is for an imported beam with irregular counts per slice.
  ! load_mode = "sample" levels it.

  n_thin = 0;  n_shown = 0
  do is = 1, min(stats%nslice, size(run%fbeam%slice))

    ! A slice the bunch never reached is outside the beam and not a thin load. A window
    ! that holds the slippage has a quarter of its slices empty by construction
    ! (docs/fel-physics.md sec-window). Listing those would bury the ones that matter.

    if (stats%charge_live(is, cvg%ir(1)) <= 0) cycle
    nbl = run%fbeam%slice(is)%n / max(1, run%bparam%beamlet_size)
    ratio_s = 0
    if (stats%f_angles_valid(is, cvg%ir(1)) == 1 .and. stats%f_p_in(is, cvg%ir(1)) > 0) &
          ratio_s = max(0.0_rp, stats%f_power(is, cvg%ir(1)) - stats%f_p_in(is, cvg%ir(1))) / &
                    stats%f_p_in(is, cvg%ir(1))
    if (nbl >= 128 .and. ratio_s <= 0.1_rp) cycle
    n_thin = n_thin + 1
    if (n_shown >= n_thin_max$) cycle
    if (n_shown == 0) call out_io (s_blank$, r_name, '             Thin slices at the mode ' // &
          'peak, under 128 beamlets or over a tenth outside the mode:')
    n_shown = n_shown + 1
    write (line, '(a, i0, a, 2a, i0, a, es9.2)') '               slice ', is, ', z = ', &
          trim(adjustl(fel_si_str((is - 1) * run%fbeam%slice_spacing, 'm'))), ', ', nbl, &
          ' beamlets, ratio ', ratio_s
    call out_io (s_blank$, r_name, trim(line))
  enddo
  if (n_thin > n_shown) then
    write (line, '(a, i0, a, i0, a, i0)') '               and ', n_thin - n_shown, ' more, ', &
          n_thin, ' thin of ', stats%nslice
    call out_io (s_blank$, r_name, trim(line))
  endif
  if (n_thin > 0) call out_io (s_blank$, r_name, '             Possible solution: ' // &
        'load_mode = "sample", which draws the same number of beamlets in every slice.')
endif

! Where the time went, one row per phase the run entered. The phases partition the
! walk. So the unaccounted row is the walk's own overhead, and nothing is counted twice.
! docs/performance.md reads a table like this one.

call fel_timer_write (r_name)

! The estimate beside the clock, on every run, so the two can be compared on any
! machine (fel_cost_mod).

call fel_cost_estimate (run, secs, backend)
write (line, '(7a)') ' Cost        walk ', &
      trim(fel_cost_secs_str(fel_timer_seconds(fel_t_walk$))), ' measured, ', &
      trim(fel_cost_secs_str(secs)), ' estimated on the ' // trim(backend) // ' path for ', &
      fel_cost_machine$, '.'
call out_io (s_blank$, r_name, trim(line))

! The FP32 lockstep instrument's worst case, when it ran. The full per-step record
! and the summary block are in the .fp32.txt file listed below.

if (run%fp32%on) then
  write (line, '(a, 3(a, es9.2), a, es9.2)') ' FP32        worst', &
        ' pz ', run%fp32%worst(5), ', theta ', run%fp32%worst(6), &
        ', phasor ', run%fp32%worst(7), ', guard ulp ', run%fp32%ulp_min
  call out_io (s_blank$, r_name, trim(line))
endif

! What slippage carried out of the window, per member of the field set. The energy the
! record lost this way is the other half of a conservation statement the window's own
! power cannot make. Otherwise the unaveraged ledger is its only reader.

if (run%ffield(1)%slip%timerun) then
  line = ' Escaped    '
  do ih = 1, run%n_harm
    ! Full precision, since a restart comparison differences this against the same total
    ! at a checkpoint (docs/validation.md). Six digits cannot reach its floor.
    write (line, '(2a, i0, a, es24.16, a)') trim(line), ' h', run%ffield(ih)%harm, ' ', &
          run%ffield(ih)%slip%u_escaped, ' J'
    if (ih < run%n_harm) line = trim(line) // ','
  enddo
  call out_io (s_blank$, r_name, trim(line) // ' transmitted out of the window by slippage')
endif

! The file list is by existence, and not by replaying which switches were on. A
! candidate name that is there is listed with its size. The naming rules stay in the
! writers, which own them.

n_listed = 0
call note_file (trim(run%global%out_root) // '-final.beam.h5')
do ih = 1, run%n_harm
  hsuf = ''
  if (run%ffield(ih)%harm /= 1) write (hsuf, '(a, i0)') '-h', run%ffield(ih)%harm
  call note_file (trim(run%global%out_root) // '-final' // trim(hsuf) // '.wf.h5')
  call note_file (trim(run%global%out_root) // '-escaped' // trim(hsuf) // '.fld.h5')
  call note_file (trim(run%global%out_root) // '-pulse' // trim(hsuf) // '.fld.h5')
enddo
call note_file (trim(run%global%out_root) // '.stats.h5')
call note_file (trim(run%global%out_root) // '.diag.txt')
call note_file (trim(run%global%out_root) // '.ledger.txt')
call note_file (trim(run%global%out_root) // '.import.txt')
call note_file (trim(run%global%out_root) // '.fp32.txt')
call note_file (trim(run%global%out_root) // '.migration.txt')

call out_io (s_blank$, r_name, &
      '================================================================================')

!------------------------------------------------------------------------------
contains

!+
! Subroutine note_file (name)
!
! Routine to list one output file with its size, if it is there. Existence is the test,
! so the footer never claims a file the run did not write.
!-

subroutine note_file (name)

character(*) name
character(300) row
integer(8) fsize
logical there

!

inquire (file = name, exist = there, size = fsize)
if (.not. there) return
n_listed = n_listed + 1
if (n_listed == 1) then
  write (row, '(2a, t62, a)') ' Wrote       ', trim(name), trim(adjustl(fel_si_str(real(fsize, rp), 'B')))
else
  write (row, '(2a, t62, a)') '             ', trim(name), trim(adjustl(fel_si_str(real(fsize, rp), 'B')))
endif
call out_io (s_blank$, r_name, trim(row))

end subroutine note_file

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
logical rerr, err_flag, has_res
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
ffield(ihh)%wf%dz = fbeam%slice_spacing
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
  if (.not. ieee_is_finite(residual) .or. (fbeam%n_wavelength >= 1 .and. abs(residual) >= fbeam%n_wavelength)) then
    call h5fclose_f (f_id, h5_err)
    call out_io (s_error$, r_name, 'THE FIELD FILE''S slippageResidual IS NOT A VALUE THIS TRACKER WRITES: \es12.4\ ', &
                 'ITS MAGNITUDE MUST STAY BELOW n_wavelength = \i0\. THE FILE: ' // trim(fname), &
                 r_array = [residual], i_array = [fbeam%n_wavelength])
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
       'BEFORE THE GROUP EXISTED IS NOT ONE, AND NOTHING IS TAKEN FROM IT (docs/reading-output.md).')
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
! Routine to drain the slippage bank. Each slice carried out streams to
! <out_root>-escaped.fld.h5 as it leaves, so peak memory is a handful of grid planes.
! Each is stored with its wavefront_params and its transmission z. The params are the
! full 4x4 at bank time, which the analytic free-space propagation of pulse statistics
! needs. Genesis field-file conventions (dfl units) are used, so existing tools read
! the file.
!
! During the run the groups are numbered in transmission order: slice k is the k-th to
! leave. That is the reverse of Genesis's order. Light leaves at the head, so the first
! to leave is the furthest ahead. fel_order_bank renumbers them at finalize, and the
! root datasets are written then. So the finished file is a Genesis field file.
!
! Input:
!   run       -- fel_run_struct: Run state. run%ffield(ihh)%bank holds the slices to drain.
!   ihh       -- integer: Field-set index whose bank is drained.
!
! Output:
!   run       -- fel_run_struct: Banked params appended to run%bank_z / run%bank_pms.
!                  run%n_banked(ihh) and run%esc_id(ihh) updated.
!   err_flag  -- logical: Set True on a write error. False otherwise.
!-

subroutine fel_drain_bank (run, ihh, err_flag)

type (fel_run_struct), target :: run
type (fel_field_struct), pointer :: ffield(:)
type (fel_beam_struct), pointer :: fbeam
real(rp), pointer :: z_now
logical two_pol, err_flag
integer n_harm

type (wavefront_params_struct) pms
type (fel_bank_struct), pointer :: bnk
type (wavefront_struct), pointer :: wfl
integer(hid_t) g_id
real(rp), allocatable :: work(:), gz(:,:), gp(:,:,:)
real(rp) dfl_scale
integer ihh, k, nx, h5e
logical berr
character(20) gname

!

err_flag = .true.
ffield => run%ffield
fbeam => run%fbeam
z_now => run%z_now
two_pol = run%two_pol
n_harm = run%n_harm

bnk => ffield(ihh)%bank
wfl => ffield(ihh)%wf
if (bnk%n == 0) then
  err_flag = .false.
  return
endif
nx = size(wfl%Ex, 1)
dfl_scale = wfl%dx / sqrt(2 * (mu_0_vac * c_light))
allocate (work(nx*nx))

if (run%esc_id(ihh) == 0) then
  call hdf5_open_file (trim(fel_escaped_file_name(run, ihh)), 'WRITE', run%esc_id(ihh), berr)
  if (berr) return
endif

do k = 1, bnk%n
  run%n_banked(ihh) = run%n_banked(ihh) + 1
  if (.not. allocated(run%bank_z)) then
    allocate (run%bank_z(1024, n_harm), run%bank_pms(25, 1024, n_harm))
  elseif (run%n_banked(ihh) > size(run%bank_z, 1)) then
    call move_alloc (run%bank_z, gz);  call move_alloc (run%bank_pms, gp)
    allocate (run%bank_z(2*size(gz,1), n_harm), run%bank_pms(25, 2*size(gz,1), n_harm))
    run%bank_z(1:size(gz,1), :) = gz;  run%bank_pms(:, 1:size(gz,1), :) = gp
    deallocate (gz, gp)
  endif

  call wavefront_params_of_plane (bnk%plane(:,:,k), wfl%dx, wfl%wavelength, &
                                  fbeam%slice_spacing, pms, .true., berr)
  if (berr) return
  pms%s = z_now
  run%bank_z(run%n_banked(ihh), ihh) = z_now
  run%bank_pms(:, run%n_banked(ihh), ihh) = [pms%centroid, reshape(pms%sigma, [16]), &
                           pms%energy, pms%power, pms%on_axis_intensity, pms%emit_x, pms%emit_y]
  if (two_pol) then             ! The y component's params, written below as wavefront_params_y.
    call wavefront_params_of_plane (bnk%plane_y(:,:,k), wfl%dx, wfl%wavelength, &
                                    fbeam%slice_spacing, pms, .true., berr)
    if (berr) return
  endif

  write (gname, '(a, i0.6)') 'slice', run%n_banked(ihh)
  call H5Gcreate_f (run%esc_id(ihh), trim(gname), g_id, h5e)
  if (h5e < 0) return
  work = dfl_scale * reshape(real(bnk%plane(:,:,k), rp), [nx*nx])
  call hdf5_write_dataset_real (g_id, 'field-real', work, berr);  if (berr) return
  work = dfl_scale * reshape(aimag(bnk%plane(:,:,k)), [nx*nx])
  call hdf5_write_dataset_real (g_id, 'field-imag', work, berr);  if (berr) return
  if (two_pol) then
    work = dfl_scale * reshape(real(bnk%plane_y(:,:,k), rp), [nx*nx])
    call hdf5_write_dataset_real (g_id, 'field-real-y', work, berr);  if (berr) return
    work = dfl_scale * reshape(aimag(bnk%plane_y(:,:,k)), [nx*nx])
    call hdf5_write_dataset_real (g_id, 'field-imag-y', work, berr);  if (berr) return
    call hdf5_write_dataset_real (g_id, 'wavefront_params_y', &
            [pms%centroid, reshape(pms%sigma, [16]), pms%energy, pms%power, &
             pms%on_axis_intensity, pms%emit_x, pms%emit_y], berr);  if (berr) return
  endif
  call hdf5_write_dataset_real (g_id, 'z_transmit', [z_now], berr);  if (berr) return
  call hdf5_write_dataset_real (g_id, 'wavefront_params', run%bank_pms(:, run%n_banked(ihh), ihh), berr);  if (berr) return
  call H5Gclose_f (g_id, h5e)
enddo

err_flag = .false.

end subroutine fel_drain_bank

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_order_bank (run, ihh, err_flag)
!
! Routine to renumber the escaped file's slice groups into Genesis's order, where
! slice i sits at refposition + (i-1) slicespacing. fel_drain_bank numbers them in
! transmission order as they leave. Light leaves at the window's head, so the k-th to
! leave is the (n_banked - k + 1)-th in increasing coordinate. Each group moves whole,
! with its field planes, wavefront_params and z_transmit. This is done once, at
! finalize, when n_banked is known. run%bank_z and run%bank_pms stay in transmission
! order.
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
integer ihh, k, n, h5e
logical err_flag
character(20) gname, tname

!

err_flag = .true.
n = run%n_banked(ihh)

! Two passes, through names no slice ever carries, since a direct move would land on a
! group that still exists.

do k = 1, n
  write (gname, '(a, i0.6)') 'slice', k
  write (tname, '(a, i0.6)') 'leaving', k
  call H5Lmove_f (run%esc_id(ihh), trim(gname), run%esc_id(ihh), trim(tname), h5e)
  if (h5e < 0) return
enddo

do k = 1, n
  write (tname, '(a, i0.6)') 'leaving', k
  write (gname, '(a, i0.6)') 'slice', n - k + 1
  call H5Lmove_f (run%esc_id(ihh), trim(tname), run%esc_id(ihh), trim(gname), h5e)
  if (h5e < 0) return
enddo

err_flag = .false.

end subroutine fel_order_bank

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_escaped_file_name (run, ihh) result (fname)
!
! Routine to construct the escaped-field file name of field ihh. The fundamental
! keeps its pre-harmonic name. A harmonic's file carries -h<h>.
!
! Input:
!   run    -- fel_run_struct: Run state (out_root and the field set).
!   ihh    -- integer: Field-set index.
!
! Output:
!   fname  -- character(420): The escaped-field file name.
!-

function fel_escaped_file_name (run, ihh) result (fname)

type (fel_run_struct), target :: run
type (fel_field_struct), pointer :: ffield(:)
character(400) out_root
integer ihh
character(420) fname

ffield => run%ffield
out_root = run%global%out_root

if (ffield(ihh)%harm == 1) then
  fname = trim(out_root) // '-escaped.fld.h5'
else
  write (fname, '(2a, i0, a)') trim(out_root), '-escaped-h', ffield(ihh)%harm, '.fld.h5'
endif

end function fel_escaped_file_name

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_finalize_diagnostics (run, err_flag)
!
! Routine to finalize the run's file diagnostics: the stats file always. With
! keep_escaped_field, also the escaped file's root datasets and the full pulse at the
! exit plane.
!
! For the pulse, each banked slice is read back and free-space propagated over
! z_end - z_transmit. Transmitted light is fixed information, and undulator vacuum is
! free space for light with no beam under it. The slices are concatenated above the
! live window. Earliest-transmitted light is furthest ahead. So banked slice k lands at
! pulse index nslice + (run%n_banked - k + 1), which is also its number in the
! finished escaped file (fel_order_bank). The caller has already unrotated the live
! window (the final-dump cshift).
!
! refposition in both files is what Genesis means by it: the in-bunch coordinate of
! the file's slice 1 (readFieldHDF5 places slice i at refposition + (i-1) slicespacing).
! It is 0 for the pulse, whose slice 1 is the live window's tail. It is nslice times the
! slice spacing for the escaped file, whose slice 1 is the slot just past the window's
! head. Neither is the beamline position, which is the openPMD files' zCoordinate.
!
! Input:
!   run       -- fel_run_struct: Run state at the end of tracking.
!
! Output:
!   run       -- fel_run_struct: Escaped files closed (run%esc_id zeroed).
!   err_flag  -- logical: Set True on a write error. False otherwise.
!-

subroutine fel_finalize_diagnostics (run, err_flag)

type (fel_run_struct), target :: run
type (fel_field_struct), pointer :: ffield(:)
type (fel_beam_struct), pointer :: fbeam
type (fel_stats_params_struct) sprm
integer nx, h5e, ihh, n_harm, nslice, is
logical ferr, err_flag, two_pol, keep_escaped_field
character(400) out_root
character(420) fname_h
character(*), parameter :: r_name = 'fel_finalize_diagnostics'

!

err_flag = .true.
ffield => run%ffield
fbeam => run%fbeam
n_harm = run%n_harm
nslice = run%nslice
two_pol = run%two_pol
keep_escaped_field = run%global%keep_escaped_field
out_root = run%global%out_root

! The parameters the file states as data. They are assembled here because the whole run
! is visible here. fel_stats_mod knows the accumulator and nothing else.

sprm%slice_spacing = fbeam%slice_spacing
sprm%species = 'electron'
sprm%beta0 = fel_p0_mc(fbeam) / sqrt(fel_p0_mc(fbeam)**2 + 1)
sprm%grid_dx = run%ffield(1)%wf%dx
sprm%grid_dy = run%ffield(1)%wf%dy

call fel_stats_write (run%stats, sprm, trim(out_root) // '.stats.h5', ferr)
call fel_write_params (run, trim(out_root) // '.stats.h5')
if (ferr) return
call fel_write_lattice (run, trim(out_root) // '.stats.h5')
call fel_write_meta (run, trim(out_root) // '.stats.h5')

if (.not. keep_escaped_field) then
  err_flag = .false.
  return
endif

do ihh = 1, n_harm
  nx = size(ffield(ihh)%wf%Ex, 1)

  if (run%esc_id(ihh) /= 0) then
    call fel_order_bank (run, ihh, ferr)
    if (ferr) return
    call hdf5_write_dataset_int  (run%esc_id(ihh), 'gridpoints',   [nx],                        ferr)
    call hdf5_write_dataset_real (run%esc_id(ihh), 'gridsize',     [ffield(ihh)%wf%dx],         ferr)
    call hdf5_write_dataset_real (run%esc_id(ihh), 'refposition',  [nslice * fbeam%slice_spacing], ferr)
    call hdf5_write_dataset_real (run%esc_id(ihh), 'wavelength',   [ffield(ihh)%wf%wavelength], ferr)
    call hdf5_write_dataset_int  (run%esc_id(ihh), 'slicecount',   [run%n_banked(ihh)],             ferr)
    call hdf5_write_dataset_real (run%esc_id(ihh), 'slicespacing', [fbeam%slice_spacing],       ferr)
    call H5Fclose_f (run%esc_id(ihh), h5e)
    run%esc_id(ihh) = 0
    call out_io (s_blank$, r_name, '  ' // trim(fel_escaped_file_name(run, ihh)))
  endif

  ! The full pulse at the exit plane. With two live polarizations, one file per
  ! component (Genesis's format holds one). Harmonic pulses carry -h<h>.

  if (two_pol) then
    call fel_write_pulse_file (run, trim(out_root) // '-pulse-x.fld.h5', .false., ihh, ferr)
    if (ferr) return
    call fel_write_pulse_file (run, trim(out_root) // '-pulse-y.fld.h5', .true., ihh, ferr)
    if (ferr) return
    call out_io (s_blank$, r_name, '  ' // trim(out_root) // '-pulse-{x,y}.fld.h5')
  elseif (ffield(ihh)%harm == 1) then
    call fel_write_pulse_file (run, trim(out_root) // '-pulse.fld.h5', .false., ihh, ferr)
    if (ferr) return
    call out_io (s_blank$, r_name, '  ' // trim(out_root) // '-pulse.fld.h5')
  else
    write (fname_h, '(2a, i0, a)') trim(out_root), '-pulse-h', ffield(ihh)%harm, '.fld.h5'
    call fel_write_pulse_file (run, trim(fname_h), .false., ihh, ferr)
    if (ferr) return
    call out_io (s_blank$, r_name, '  ' // trim(fname_h))
  endif
enddo

err_flag = .false.

end subroutine fel_finalize_diagnostics

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_write_pulse_file (run, fname, use_y, ihh, err_flag)
!
! Routine to write one component of the full exit-plane pulse (Genesis's format holds
! one per file). It writes the live window's slices, then each banked slice read back
! from the escaped file and free-space propagated over z_end - z_transmit. Shared by
! the -x and -y files. The escaped file is in Genesis's order by now (fel_order_bank).
! So its slice i lands at pulse slice nslice + i. The transmission z of that slice is
! run%bank_z(n_banked - i + 1), which stays in transmission order.
!
! Input:
!   run       -- fel_run_struct: Run state.
!   fname     -- character(*): Pulse file to write.
!   use_y     -- logical: If True write the Ey component, otherwise the Ex component.
!   ihh       -- integer: Field-set index.
!
! Output:
!   err_flag  -- logical: Set True on a read or write error. False otherwise.
!-

subroutine fel_write_pulse_file (run, fname, use_y, ihh, err_flag)

type (fel_run_struct), target :: run
type (fel_field_struct), pointer :: ffield(:)
type (fel_beam_struct), pointer :: fbeam
type (wavefront_struct) wf1
real(rp), pointer :: z_now
character(*) fname
logical use_y, err_flag
integer ihh, n_harm, nslice
type (wavefront_struct), pointer :: wfl
integer(hid_t) p_id, e_id, g_id
real(rp), allocatable :: work(:), re_w(:), im_w(:)
real(rp) dfl_scale
integer k, is_f, nx, h5e
logical ferr
character(24) gname, dset_r, dset_i

!

err_flag = .true.
ffield => run%ffield
fbeam => run%fbeam
z_now => run%z_now
n_harm = run%n_harm
nslice = run%nslice

wfl => ffield(ihh)%wf
nx = size(wfl%Ex, 1)
dfl_scale = wfl%dx / sqrt(2 * (mu_0_vac * c_light))
dset_r = 'field-real';  dset_i = 'field-imag'
if (use_y) then
  dset_r = 'field-real-y';  dset_i = 'field-imag-y'
endif

call hdf5_open_file (trim(fname), 'WRITE', p_id, ferr);  if (ferr) return
call hdf5_write_dataset_int  (p_id, 'gridpoints',   [nx],                       ferr)
call hdf5_write_dataset_real (p_id, 'gridsize',     [wfl%dx],                   ferr)
call hdf5_write_dataset_real (p_id, 'refposition',  [0.0_rp],                   ferr)
call hdf5_write_dataset_real (p_id, 'wavelength',   [wfl%wavelength],           ferr)
call hdf5_write_dataset_int  (p_id, 'slicecount',   [nslice + run%n_banked(ihh)],   ferr)
call hdf5_write_dataset_real (p_id, 'slicespacing', [fbeam%slice_spacing],      ferr)

allocate (work(nx*nx), re_w(nx*nx), im_w(nx*nx))

do is_f = 1, nslice
  write (gname, '(a, i0.6)') 'slice', is_f
  call H5Gcreate_f (p_id, trim(gname), g_id, h5e)
  if (use_y) then
    work = dfl_scale * reshape(real(wfl%Ey(:,:,is_f), rp), [nx*nx])
  else
    work = dfl_scale * reshape(real(wfl%Ex(:,:,is_f), rp), [nx*nx])
  endif
  call hdf5_write_dataset_real (g_id, 'field-real', work, ferr);  if (ferr) return
  if (use_y) then
    work = dfl_scale * reshape(aimag(wfl%Ey(:,:,is_f)), [nx*nx])
  else
    work = dfl_scale * reshape(aimag(wfl%Ex(:,:,is_f)), [nx*nx])
  endif
  call hdf5_write_dataset_real (g_id, 'field-imag', work, ferr);  if (ferr) return
  call H5Gclose_f (g_id, h5e)
enddo

if (run%n_banked(ihh) > 0) then
  if (.not. allocated(wf1%Ex)) then
    allocate (wf1%Ex(nx, nx, 1))
    wf1%dx = wfl%dx;  wf1%dy = wfl%dy;  wf1%dz = fbeam%slice_spacing
  endif
  wf1%wavelength = wfl%wavelength      ! Per field: banked light drifts at its wavelength.

  call hdf5_open_file (trim(fel_escaped_file_name(run, ihh)), 'READ', e_id, ferr);  if (ferr) return
  do k = 1, run%n_banked(ihh)
    write (gname, '(a, i0.6)') 'slice', k
    g_id = hdf5_open_group (e_id, trim(gname), ferr, .true.);  if (ferr) return
    call hdf5_read_dataset_real (g_id, trim(dset_r), re_w, ferr, trim(gname));  if (ferr) return
    call hdf5_read_dataset_real (g_id, trim(dset_i), im_w, ferr, trim(gname));  if (ferr) return
    call H5Gclose_f (g_id, h5e)
    wf1%Ex(:,:,1) = reshape(cmplx(re_w, im_w, wf_rp), [nx, nx]) / dfl_scale

    call wavefront_drift (wf1, z_now - run%bank_z(run%n_banked(ihh) - k + 1, ihh), ferr);  if (ferr) return

    write (gname, '(a, i0.6)') 'slice', nslice + k
    call H5Gcreate_f (p_id, trim(gname), g_id, h5e)
    work = dfl_scale * reshape(real(wf1%Ex(:,:,1), rp), [nx*nx])
    call hdf5_write_dataset_real (g_id, 'field-real', work, ferr);  if (ferr) return
    work = dfl_scale * reshape(aimag(wf1%Ex(:,:,1)), [nx*nx])
    call hdf5_write_dataset_real (g_id, 'field-imag', work, ferr);  if (ferr) return
    call H5Gclose_f (g_id, h5e)
  enddo
  call H5Fclose_f (e_id, h5e)
endif

call H5Fclose_f (p_id, h5e)

err_flag = .false.

end subroutine fel_write_pulse_file

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
! and text, how many files the parser opened, an ISO timestamp and the Bmad version.
!
! The texts are datasets and not attributes. HDF5 caps a single attribute at 64 kB, and
! the largest that writes is 65495 bytes. A scalar string dataset takes 3 MB. An echoed
! deck runs to 12 kB and a lattice text to 37 kB on a real lattice
! (docs/reading-output.md). A failure would only warn, so the provenance would go missing
! silently on a lattice of ordinary size. Since input_echo carries out_root, two runs
! differ in meta/. The benchmark excludes meta/ by name from identity comparisons.
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
! timestamp and the Bmad version identify the run. The login name and working directory
! go in only under run.record_environment. The file records a lattice basename.
! The deck is echoed as read, so an absolute path there is still an absolute path
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
logical merr
character(16) table
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
      'What the user controlled: one subgroup per table of the deck, every key ' // &
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
    call H5LTset_attribute_string_f (b_id, 'resample', 'description', 'The deck''s [beam.resample] table.', h5e)
  else
    call H5Gcreate_f (p_id, trim(table), g_id, h5e)
    call H5LTset_attribute_string_f (p_id, trim(table), 'kind', 'input', h5e)
    call H5LTset_attribute_string_f (p_id, trim(table), 'description', &
          'The deck''s [' // trim(table) // '] table.', h5e)
  endif
  merr = merr .or. (h5e < 0)
  do i = 1, size(keys)
    if (keys(i)%table /= table .or. .not. keys(i)%recorded) cycle
    call write_key (keys(i))
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
        ' As the deck wrote it.', sv, merr)
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
! Routine to render the resolved deck as one string, every key of the schema with its
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
