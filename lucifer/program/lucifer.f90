!+
! Program fel_track_test
!
! FEL tracker validated against Genesis 1.3 Version 4: single-slice steady state
! (deliverable 3) and multi-slice time dependence with slippage (deliverable 4). See
! bsim/fel/README.md and the design brief.
!
! The program walks a Bmad lattice and applies the seam of the design (brief section 4.1):
!
!   - Elements whose name starts with UND are FEL segments. They are stepped internally in
!     delz with the transcribed Genesis physics (fel_track_mod): transverse push with
!     natural focusing, RK4 ponderomotive advance, source deposition and FFT field solve.
!     Bmad tracking is not used inside them.
!
!   - Every other element: each slice's bunch is converted from the packed FEL arrays to
!     coord_structs and tracked by Bmad (track1_bunch). The packed arrays ARE Bmad
!     coordinates (see fel_beam_mod), so the conversion is a plain copy; the only phase
!     bookkeeping is one advance of the common reference phase phi0 per element. The
!     radiation field is drifted through free space by wavefront_drift.
!
! Time dependence follows from the starting dumps alone: a multi-slice dump makes a
! time-dependent run (slippage active), a single-slice dump the steady state, with no
! separate switch -- the same rule as Genesis, whose imports carry the time window
! (ImportBeam.cpp). Slippage is precomputed as a schedule over the lattice, transcribing
! Lattice::calcSlippage with Genesis's reference gamma:
!
!   - each undulator step slips dz*(1+aw^2)/(2*gamma0^2*lambda) wavelengths;
!   - interlude elements slip zero, but their lengths accumulate, and when an undulator
!     follows, floor(Lz/(2*gamma0^2*lambda)) + 1 wavelengths of autophasing land on the
!     last interlude element before it (Lattice.cpp:173: "auto phasing would always add
!     some slippage" -- the field record shifts an integer number of wavelengths, the
!     particle phase does not move, exactly as in Genesis where the fractional part is
!     commented out);
!   - a trailing interlude section adds the same fixup to the lattice's last element
!     (Lattice.cpp:193).
!
! The schedule is applied through fel_apply_slippage after each step's field solve and
! before its diagnostics -- Gencore's step order (track beam, track field, slippage,
! diagnostics). The field record rotates rather than moves; everything reading it in time
! order (the per-slice diagnostics here, the final dump) goes through fel_field_index.
!
! The starting state is a pair of Genesis dumps (&write of beam and field), so both codes
! track from bitwise-identical initial conditions. Diagnostics matching Genesis's
! definitions are recorded at the same z positions Genesis records them: once at the start
! and once after every integration step, one step per interlude element -- one row per
! slice per record, in time-window order.
!
! Input is a namelist file:
!
!   &fel_track_params
!     lat_file = "aramis.bmad"                 ! Bmad lattice.
!     beam_file = "Aramis-initial.par.h5"      ! Genesis particle dump to start from.
!     field_file = "Aramis-initial.fld.h5"     ! Genesis field dump to start from.
!     out_root = "fel_td"                      ! Prefix for the three output files.
!     gamma0 = 11357.82                        ! Genesis's reference gamma.
!     delz = 0.045                             ! Target integration step inside undulators [m].
!     und_aw = 0.84853                         ! Undulator parameter (rms).
!     und_lambdau = 0.015                      ! Undulator period [m].
!     und_kx = 0.5, und_ky = 0.5               ! Natural focusing, deck convention (before ku^2).
!     und_helical = T
!     interlude_model = "bmad"                 ! "bmad" (the seam, default) or "genesis".
!     split_weights = F                        ! Weight-invariance test mode; see below.
!   &end
!
! interlude_model selects how the field-free elements are handled. "bmad" is the
! deliverable's architecture: track1_bunch for the particles, the exact theta mapping from
! Bmad's z, wavefront_drift for the field. "genesis" instead uses the transcribed Genesis
! interlude step (fel_track_interlude_genesis) everywhere, which prices what the seam
! changes: with it the whole run should agree with Genesis at transcription level, and the
! difference between the two modes is the transport model difference, measured rather than
! argued about. The slippage schedule is identical in both models.
!
! split_weights = T replaces each imported particle by two copies at identical
! coordinates carrying 1/3 and 2/3 of its weight. Every collective observable -- power,
! bunching, the field itself -- must be identical to the unsplit run, because the
! dynamics is per particle and the sources and reductions are linear in the weight. The
! benchmark harness runs this against the unsplit run to test the weighted paths, which
! nothing Genesis produces can test: the Genesis dump format carries no weights, so a
! Genesis comparison only ever sees the uniform case.
!
! Outputs: <out_root>.diag.txt (one row per slice per record: z, slice, field and beam
! diagnostics), <out_root>-final.fld.h5 and <out_root>-final.par.h5 (Genesis-format dumps
! of the end state, for field-by-field comparison; the field dump is unrotated to time
! order first, as writeFieldHDF5 does).
!-

program fel_track_test

use fel_track_mod
use wavefront_hdf5_mod
use beam_mod

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (fel_beam_struct), target :: fbeam
type (wavefront_struct) wf
type (bunch_struct) bunch
type (fel_und_struct) und
type (fel_slip_struct) slip
type (fel_slice_diag_struct) bdiag

real(rp) :: gamma0 = 0, delz = 0, und_aw = 0, und_lambdau = 0
real(rp) :: und_kx = 0.5_rp, und_ky = 0.5_rp
logical :: und_helical = .true.
logical :: split_weights = .false.
character(400) :: lat_file = '', beam_file = '', field_file = '', out_root = 'fel_track'
character(16) :: interlude_model = 'bmad'

real(rp), allocatable :: ele_slip(:)     ! Slippage applied after each element's last step [wavelengths].
real(rp) z_now, ks, qf, und_slip_step, Lz, gamma0_ref
integer ie, is, istep, n_arg, iu_diag, iu_nml, nslice, prev_ie
logical err, timerun

character(400) param_file
character(*), parameter :: r_name = 'fel_track_test'

namelist / fel_track_params / lat_file, beam_file, field_file, out_root, gamma0, delz, &
                           und_aw, und_lambdau, und_kx, und_ky, und_helical, interlude_model, split_weights

! Read parameters.

n_arg = command_argument_count()
if (n_arg /= 1) then
  print '(a)', 'Usage: fel_track_test <param_file>'
  stop 1
endif
call get_command_argument (1, param_file)

open (newunit = iu_nml, file = param_file, status = 'old', action = 'read')
read (iu_nml, nml = fel_track_params)
close (iu_nml)

if (gamma0 <= 0 .or. delz <= 0 .or. und_aw <= 0 .or. und_lambdau <= 0) then
  print '(a)', 'fel_track_test: gamma0, delz, und_aw and und_lambdau must all be set and positive.'
  stop 1
endif

if (interlude_model /= 'bmad' .and. interlude_model /= 'genesis') then
  print '(a)', 'fel_track_test: interlude_model must be "bmad" or "genesis", got: ' // trim(interlude_model)
  stop 1
endif

! Read the lattice and the shared starting state.

call bmad_parser (lat_file, lat)
branch => lat%branch(0)

call fel_read_genesis4_beam (fbeam, beam_file, gamma0, err)
if (err) stop 1
if (split_weights) call do_split_weights (fbeam)
nslice = size(fbeam%slice)

call wavefront_read_genesis4 (wf, field_file, err)
if (err) stop 1
ks = twopi / wf%wavelength

! The beam and field dumps must describe the same time window: one field slice per beam
! slice, at the same wavelength. Checked, never assumed (FINDINGS.md section 5).

if (size(wf%Ex, 3) /= nslice) then
  print '(2(a, i0))', 'fel_track_test: beam has ', nslice, ' slices but the field has ', size(wf%Ex, 3)
  stop 1
endif
if (abs(wf%wavelength - fbeam%wavelength) > 1e-12_rp * fbeam%wavelength) then
  print '(a, 2es20.12)', 'fel_track_test: beam and field dumps disagree on the wavelength: ', &
                         fbeam%wavelength, wf%wavelength
  stop 1
endif

! Time dependence follows from the dumps: more than one slice makes a time-dependent run
! with slippage active; one slice is the steady state and fel_apply_slippage is a no-op.

timerun = (nslice > 1)
slip%timerun = timerun
slip%sample = fbeam%slice_spacing / fbeam%wavelength
gamma0_ref = fel_gamma0(fbeam)

! Undulator segment parameters, constant for every segment in this benchmark. kx, ky get
! Genesis's unroll scaling by ku^2 (Lattice.cpp:412-413).

und%aw = und_aw
und%ku = twopi / und_lambdau
und%kx = und_kx * und%ku**2
und%ky = und_ky * und%ku**2
und%helical = und_helical

! Undulator slippage per integration step, in wavelengths: dz/(2*gamma0^2*lambda/(1+aw^2))
! (Lattice.cpp:167-168 with Genesis's reference gamma). Note the step length is set per
! element below; this is the rate per meter, multiplied by dz at use.

und_slip_step = (1 + und_aw**2) / (2 * gamma0_ref**2 * wf%wavelength)

! The rest of the schedule: drift autophasing. Interludes accumulate Lz; the last
! interlude before each undulator gets floor(Lz/(2*gamma0^2*lambda)) + 1 wavelengths
! (Lattice.cpp:171-174, guarded there by Lz > 0). The end-of-lattice fixup
! (Lattice.cpp:191-193) is UNGUARDED in Genesis: the last element always gets
! floor(Lz/(2*gamma0^2*lambda)) + 1, which is +1 even with no trailing interlude at all
! ("autophasing is applied in case for [a] second, succeeding run"). Transcribed as is --
! omitting that +1 leaves the field record one rotation short at the very end, found the
! hard way against the single-segment time-dependent run.

allocate (ele_slip(branch%n_ele_track))
ele_slip = 0
Lz = 0
prev_ie = 0

do ie = 1, branch%n_ele_track
  ele => branch%ele(ie)
  if (ele%value(l$) == 0) cycle
  if (ele%name(1:3) == 'UND') then
    if (Lz > 0 .and. prev_ie > 0) then
      ele_slip(prev_ie) = ele_slip(prev_ie) + floor(Lz / (2 * gamma0_ref**2 * wf%wavelength)) + 1
      Lz = 0
    endif
  else
    Lz = Lz + ele%value(l$)
  endif
  prev_ie = ie
enddo
if (prev_ie > 0) then
  ele_slip(prev_ie) = ele_slip(prev_ie) + floor(Lz / (2 * gamma0_ref**2 * wf%wavelength)) + 1
endif

! Diagnostics file, one row per slice per record at Genesis's record positions, slices in
! time-window order.

open (newunit = iu_diag, file = trim(out_root) // '.diag.txt', action = 'write')
write (iu_diag, '(a, i0)') '# nslice = ', nslice
write (iu_diag, '(a)') '#         z            slice        power         on_axis_intensity        bunching        ' // &
      'bunching_phase        mean_gamma          sigma_gamma           sigma_x               sigma_y'

z_now = 0
call write_diag_rows()     ! Initial record, matching Genesis's diag before the first step.

! Walk the lattice.

do ie = 1, branch%n_ele_track
  ele => branch%ele(ie)

  ! Zero length elements (Bmad's end marker, for one) get no step and no diagnostic
  ! record: Genesis's unrolled lattice has no counterpart for them.

  if (ele%value(l$) == 0) cycle

  if (ele%name(1:3) == 'UND') then

    ! FEL segment: Genesis's unroll, nstep = round(l/delz), equal steps. Slippage after
    ! each step's field solve; any end-of-lattice fixup lands on the last step.

    und%nstep = nint(ele%value(l$) / delz)
    if (und%nstep == 0) und%nstep = 1
    und%dz = ele%value(l$) / und%nstep

    do istep = 1, und%nstep
      call fel_track_und_step (und, fbeam, wf, slip, err)
      if (err) stop 1
      if (istep == und%nstep) then
        call fel_apply_slippage (slip, wf, und%dz * und_slip_step + ele_slip(ie))
      else
        call fel_apply_slippage (slip, wf, und%dz * und_slip_step)
      endif
      z_now = z_now + und%dz
      call write_diag_rows()
    enddo

  elseif (interlude_model == 'bmad') then

    ! The seam: Bmad tracks each slice's bunch (coordinate copies in and out),
    ! wavefront_drift moves the field (every slice, rotation-invariant), and the common
    ! phase phi0 advances by the reference rate with Genesis's drift surrogate
    ! ks/(2*gamma0^2) as the reference wavenumber.

    do is = 1, nslice
      call fel_slice_to_bunch (fbeam, fbeam%slice(is), ele, bunch, err)
      if (err) stop 1
      call track1_bunch (bunch, ele, err)
      if (err) then
        print '(2a)', 'fel_track_test: tracking error in element ', trim(ele%name)
        stop 1
      endif
      call fel_bunch_to_slice (bunch, ele, fbeam%slice(is), err)
      if (err) stop 1
    enddo

    fbeam%phi0 = fbeam%phi0 + ele%value(l$) * &
                    fel_phi0_rate(ks, ks * 0.5_rp / gamma0_ref**2, fel_p0_mc(fbeam))

    call wavefront_drift (wf, ele%value(l$), err)
    if (err) stop 1

    call fel_apply_slippage (slip, wf, ele_slip(ie))

    z_now = z_now + ele%value(l$)
    call write_diag_rows()

  else

    ! Genesis's own interlude model, transcribed, for pricing what the seam changes.

    qf = 0
    if (ele%key == quadrupole$) qf = ele%value(k1$)
    call fel_track_interlude_genesis (qf, ele%value(l$), fbeam, wf, slip, err)
    if (err) stop 1

    call fel_apply_slippage (slip, wf, ele_slip(ie))

    z_now = z_now + ele%value(l$)
    call write_diag_rows()
  endif
enddo

close (iu_diag)

! Final dumps in Genesis format. The field record is unrotated to time order first --
! time window position is holds record slice 1 + mod(is-1+first, nslice) -- which is what
! writeFieldHDF5.cpp:86 does on the fly.

if (slip%first /= 0) then
  wf%Ex = cshift(wf%Ex, shift = slip%first, dim = 3)
  slip%first = 0
endif

call wavefront_write_genesis4 (wf, trim(out_root) // '-final.fld.h5', err, 'x')
if (err) stop 1

call fel_write_genesis4_beam (fbeam, trim(out_root) // '-final.par.h5', err)
if (err) stop 1

print '(a)', 'fel_track_test done.'
print '(a)', '  ' // trim(out_root) // '.diag.txt'
print '(a)', '  ' // trim(out_root) // '-final.fld.h5'
print '(a)', '  ' // trim(out_root) // '-final.par.h5'

call wavefront_fft_free()

!------------------------------------------------------------------------------
contains

subroutine do_split_weights (beam)

! Replace each particle by two coincident copies with weights w/3 and 2w/3. The order --
! all first copies, then all second copies -- keeps the original particles' storage
! order, which keeps the RK4 arithmetic per copy identical to the unsplit run.

type (fel_beam_struct), target :: beam
type (fel_slice_struct), pointer :: sp
integer is, ip, n0

do is = 1, size(beam%slice)
  sp => beam%slice(is)
  n0 = sp%n
  call fel_slice_reallocate (sp, 2*n0)
  do ip = 1, n0
    sp%x(n0+ip) = sp%x(ip);  sp%px(n0+ip) = sp%px(ip)
    sp%y(n0+ip) = sp%y(ip);  sp%py(n0+ip) = sp%py(ip)
    sp%z(n0+ip) = sp%z(ip);  sp%pz(n0+ip) = sp%pz(ip)
    sp%weight(n0+ip) = 2 * sp%weight(ip) / 3
    sp%weight(ip) = sp%weight(ip) / 3
  enddo
  sp%n = 2*n0
enddo

end subroutine do_split_weights

!------------------------------------------------------------------------------

subroutine write_diag_rows ()

! One row per slice, slices in time-window order: beam slice is against field slice
! fel_field_index(slip, is, nslice), the rotation Genesis applies at Field.cpp:329.

real(rp) power, on_axis
integer is

do is = 1, nslice
  call fel_field_diag (wf, fel_field_index(slip, is, nslice), power, on_axis)
  call fel_slice_diag (fbeam, fbeam%slice(is), ks, bdiag)

  write (iu_diag, '(es24.16, i8, 8es24.16)') z_now, is, power, on_axis, bdiag%bunching, &
        bdiag%bunching_phase, bdiag%mean_gamma, bdiag%sigma_gamma, bdiag%sigma_x, bdiag%sigma_y
enddo

end subroutine write_diag_rows

end program fel_track_test
