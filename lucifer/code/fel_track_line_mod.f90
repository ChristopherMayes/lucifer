!+
! Module fel_track_line_mod
!
! The walk: track_fel_line(run) is the main loop of the FEL tracker, moved verbatim
! from the driver. FEL segments are stepped with the transcribed Genesis physics,
! every other element goes through the Bmad seam (or the transcribed Genesis
! interlude model), slippage and phasing follow the precomputed schedule, and
! stats/diag records are taken at Genesis's record positions.
!
! The routine is callable repeatedly in one process (the re-entrancy contract: all
! state lives in run; the RNG is seeded at init). Errors return through err_flag;
! nothing here stops.
!-

module fel_track_line_mod

use fel_struct
use fel_io_mod
use beam_mod
use wake_mod

implicit none

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine track_fel_line (run, err_flag)
!
! Routine to walk the lattice from run%i_start through run%i_end (the resolved
! tracking window; the schedule was built on the full lattice, so a windowed run
! composes exactly with the full one). Opens and closes the run's stream files
! (diag/ledger/wake).
!
! Input:
!   run       -- fel_run_struct: Run state, fully initialized.
!
! Output:
!   run       -- fel_run_struct: Beam, fields, stats, and counters advanced to run%i_end.
!   err_flag  -- logical: Set True on any tracking or I/O error. False otherwise.
!-

subroutine track_fel_line (run, err_flag)

type (fel_run_struct), target :: run
logical err_flag

type (lat_struct), pointer :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele, wake_src
type (fel_beam_struct), pointer :: fbeam
type (fel_field_struct), pointer :: ffield(:)
type (wavefront_struct), pointer :: wf
type (fel_slip_struct), pointer :: slip
type (fel_bank_struct), pointer :: bank
type (fel_collective_struct), pointer :: coll
type (fel_stats_struct), pointer :: stats
type (fel_und_struct), pointer :: und_of(:)
integer, pointer :: fel_mode(:), fel_spp(:)
real(rp), pointer :: fel_ramp(:), ele_slip(:), fel_zoff(:), light_corr(:)
logical, pointer :: is_fel(:), dump_beam_here(:), dump_field_here(:)
type (fel_slice_diag_struct), pointer :: bdiag_arr(:)
real(rp), pointer :: fpow_arr(:), fonax_arr(:)
real(rp), pointer :: z_now, u_spont_cum, e_rad_cum, charge_dropped_tot, b_dev_max
integer, pointer :: n_moved_tot, iu_diag, iu_ledger, iu_wake

type (bunch_struct) bunch
type (fel_und_struct) und
real(rp), allocatable :: e_rad_slice(:), rad_kick(:,:)
real(rp) gamma0_ref, phase_rate, ks, qf, und_slip_step, lat_length, dE_step, dU_step, comb
integer nslice, n_harm, ie, is, ih, istep
integer(8) prog_count0, prog_count_last, prog_rate
logical write_diag, keep_escaped_field, migrate, migrate_check, any_unavg, two_pol, err
character(400) out_root
character(16) interlude_model
character(*), parameter :: r_name = 'track_fel_line'

!

err_flag = .false.
lat => run%lat
branch => lat%branch(0)
fbeam => run%fbeam
ffield => run%ffield
wf => run%ffield(1)%wf
slip => run%ffield(1)%slip
bank => run%ffield(1)%bank
coll => run%coll
stats => run%stats
und_of => run%und_of
fel_mode => run%fel_mode
fel_spp => run%fel_spp
fel_ramp => run%fel_ramp
ele_slip => run%ele_slip
fel_zoff => run%fel_zoff
light_corr => run%light_corr
is_fel => run%is_fel
dump_beam_here => run%dump_beam_here
dump_field_here => run%dump_field_here
bdiag_arr => run%bdiag_arr
fpow_arr => run%fpow_arr
fonax_arr => run%fonax_arr
z_now => run%z_now
u_spont_cum => run%u_spont_cum
e_rad_cum => run%e_rad_cum
charge_dropped_tot => run%charge_dropped_tot
b_dev_max => run%b_dev_max
n_moved_tot => run%n_moved_tot
iu_diag => run%iu_diag
iu_ledger => run%iu_ledger
iu_wake => run%iu_wake
gamma0_ref = run%gamma0_ref
phase_rate = run%phase_rate
ks = run%ks
nslice = run%nslice
n_harm = run%n_harm
any_unavg = run%any_unavg
two_pol = run%two_pol
write_diag = run%global%write_diag
keep_escaped_field = run%global%keep_escaped_field
migrate = run%global%migrate
migrate_check = run%global%migrate_check
out_root = run%global%out_root
interlude_model = run%global%interlude_model
comb = run%global%comb_ds_save
run%z_last_rec = -1e30_rp

! Diagnostics file, one row per slice per record at Genesis's record positions, slices in
! time-window order.

if (write_diag) then
  open (newunit = iu_diag, file = trim(out_root) // '.diag.txt', action = 'write')
  write (iu_diag, '(a, i0)') '# nslice = ', nslice
  write (iu_diag, '(a, es22.14)') '# slice_spacing = ', fbeam%slice_spacing
  write (iu_diag, '(a)') '#         z            slice        power         on_axis_intensity        bunching        ' // &
        'bunching_phase        mean_energy         sigma_energy          sigma_x               sigma_y' // &
        '               current               n_eff'
endif

! The unaveraged energy ledger (fel-physics.tex sec:unaveraged): one row per record
! step inside FEL segments -- total weighted beam energy, total window field energy,
! and the kick-side energy change of the step. The ledger check holds
! d(E_beam + U_field) to its measured floor.

if (any_unavg) then
  open (newunit = iu_ledger, file = trim(out_root) // '.ledger.txt', action = 'write')
  write (iu_ledger, '(a)') '#         z          E_beam_rel [J]          U_field [J]           dE_kick [J]          U_escaped [J]          U_spont [J]         E_radiated [J]'
endif

n_moved_tot = 0
charge_dropped_tot = 0
b_dev_max = 0

! Progress goes to stdout, throttled by wall clock (slow modes print a steady trickle,
! fast runs just the element boundaries); scripts that redirect stdout get it as their
! log. The numbers come from the stats row just taken -- no extra computation.

call system_clock (prog_count0, prog_rate)
prog_count_last = prog_count0
lat_length = branch%ele(branch%n_ele_track)%s

z_now = branch%ele(run%i_start - 1)%s      ! 0 for a full run; the window's entry face.
if (fel_comb_take(comb, z_now, run%z_last_rec, .false.)) then
  call take_stats_record (.true.)   ! Evaluates the diag instrument too; the writer prints.
  if (err_flag) return
  call write_diag_rows()            ! Initial record, matching Genesis's diag before the first step.
endif



do ie = run%i_start, run%i_end
  ele => branch%ele(ie)

  ! Zero length elements (Bmad's end marker, for one) get no step and no diagnostic
  ! record: Genesis's unrolled lattice has no counterpart for them -- UNLESS a wake
  ! applies there (a zero-length wake element is a standard Bmad idiom), in which case
  ! the element takes the interlude wake path below. wake_src resolves the wake through
  ! lords: a wake on a superimposed or split element lives on the LORD, and ele%wake is
  ! null on its slaves (pointer_to_wake_ele; the resolution also picks exactly ONE
  ! slave of a split lord -- the one containing the lord's midpoint -- so a split wake
  ! applies once, Bmad's own convention). Checking ele%wake directly was once a real
  ! hole: lord wakes fell through to the per-slice path, where Bmad applied them
  ! within single slices and noted every zero-charge filler.

  ! With bmad_com%sr_wakes_on off, track1_bunch applies no wake and the whole-window
  ! concatenation below buys nothing but its serial cost: a wake-carrying element then
  ! routes like any other interlude (the slice loop, parallel). Measured on a user's
  ! 103-element line with 74 wake pipes at 12 threads: the switch was the difference
  ! between the serial and parallel interlude paths for the entire lattice.

  wake_src => null()
  if (bmad_com%sr_wakes_on) wake_src => pointer_to_wake_ele(ele)
  if (ele%value(l$) == 0 .and. .not. associated(wake_src)) cycle

  if (is_fel(ie)) then

    ! FEL segment: Genesis's unroll in the element's OWN step -- Bmad's standard
    ! ds_step/num_steps attributes, whose bookkeeper computes exactly Genesis's
    ! num_steps = round(l/ds_step) (attribute_bookkeeper.f90; there is no namelist
    ! step size, the same rule as every other parameter). Equal steps; slippage after
    ! each step's field solve; any end-of-lattice fixup lands on the last step.

    und = und_of(ie)
    und%bmad_transport = (fel_mode(ie) == fel_averaged$)   ! The default: bmad_standard's kernel.
    und_slip_step = (1 + und%aw**2) / (2 * gamma0_ref**2 * wf%wavelength)
    und%nstep = max(1, nint(ele%value(num_steps$)))
    und%dz = ele%value(l$) / und%nstep

    if (fel_mode(ie) == fel_unaveraged$) then
      ! The concatenated wake kick would meet the quiver-carrying chart mid-segment;
      ! nothing in this mode needs element wakes, so refuse rather than approximate.
      if (associated(wake_src)) then
        call out_io (s_error$, r_name, 'ELEMENT SR WAKES ARE NOT SUPPORTED IN THE UNAVERAGED MODE.', &
                                       'AT ELEMENT: ' // trim(ele%name))
        err_flag = .true.;  return
      endif
      call fel_unavg_setup (und, run%ustate, ele%value(l$), und%dz, fel_spp(ie), fel_ramp(ie), err)
      if (err) then
        err_flag = .true.;  return
      endif
    endif

    ! The off-phase knob (manual sec:phasing): a displaced element sees the extra
    ! upstream-break phase at entry and gives it back at exit (the downstream break
    ! is shorter by the same delta), so the anchor stays nominal downstream. Sign:
    ! positive z_offset = a longer upstream break = MORE beam delay = theta backwards,
    ! Genesis's phase-shifter convention -- anchored by the cross-code phi scan.

    if (fel_zoff(ie) /= 0) fbeam%phi0 = fbeam%phi0 - phase_rate * fel_zoff(ie)

    do istep = 1, und%nstep
      if (fel_mode(ie) == fel_unaveraged$) then
        call fel_unavg_step (und, run%ustate, fbeam, wf, slip, und%dz, istep == 1, &
                             istep == und%nstep, dE_step, dU_step, err)
        u_spont_cum = u_spont_cum + dU_step
      else
        call fel_track_und_step (und, fbeam, ffield, coll, err)
      endif
      if (err) then
        err_flag = .true.;  return
      endif
      call apply_radiation ()

      ! Element sr wake, Bmad's once-per-passage convention mirrored: one kick at the
      ! step nearest mid-element, scaled to the full element length (scale_with_length
      ! uses ele's l), applied across the WHOLE window (manual sec:seamwake). Direct
      ! kick, no transport -- this walk owns transport inside wigglers.

      if (associated(wake_src) .and. istep == (und%nstep + 1)/2) then
        call apply_bmad_wake_kick (wake_src)
        if (err_flag) return
      endif

      z_now = z_now + und%dz
      if (istep == und%nstep) then
        call apply_slippage_banked (und%dz * und_slip_step + ele_slip(ie))
        if (err_flag) return
      else
        call apply_slippage_banked (und%dz * und_slip_step)
        if (err_flag) return
      endif
      if (istep == und%nstep) call do_migrate ()
      if (err_flag) return
      if (fel_comb_take(comb, z_now, run%z_last_rec, istep == und%nstep)) then
        call take_stats_record (istep == und%nstep)
        if (err_flag) return
        if (fel_mode(ie) == fel_unaveraged$) call write_ledger_row ()
        call write_diag_rows()
      endif
      call progress_line (istep == und%nstep, istep, und%nstep)
    enddo
    if (fel_zoff(ie) /= 0) fbeam%phi0 = fbeam%phi0 + phase_rate * fel_zoff(ie)
    call end_of_element ()
    if (err_flag) return

  elseif (interlude_model == 'bmad') then

    ! The seam: Bmad tracks each slice's bunch (coordinate copies in and out),
    ! wavefront_drift moves the field (every slice, rotation-invariant), and the common
    ! phase phi0 advances by the reference rate with Genesis's drift surrogate
    ! ks/(2*gamma0^2) as the reference wavenumber.
    !
    ! This slice loop is deliberately SERIAL: track1_bunch parallelizes over particles
    ! internally (track1_bunch_hom, "$OMP parallel do if (thread_safe)", on by default
    ! via global_com%mp_threading_is_safe), so the threads are already busy inside each
    ! call, and parallelizing here as well would nest. The FEL step's parallelism over
    ! slices lives in fel_track_mod.

    if (associated(wake_src)) then

      ! Wake-carrying interlude: ALL slices as one bunch in global window coordinates,
      ! through Bmad's own track1_bunch, which applies the sr wake at ds_wake with the
      ! whole window visible head to tail (manual sec:seamwake). The per-slice path below
      ! is untouched for everything else, keeping its numerics bit-identical.

      call fel_concat_slices (fbeam, ele, run%wake_bunch, run%wake_beta0, err)
      if (err) then
        err_flag = .true.;  return
      endif
      call track1_bunch (run%wake_bunch, ele, err)
      if (err) then
        call out_io (s_error$, r_name, 'TRACKING ERROR IN ELEMENT: ' // trim(ele%name))
        err_flag = .true.;  return
      endif
      call fel_split_slices (run%wake_bunch, ele, fbeam, run%wake_beta0, .false., err)
      if (err) then
        err_flag = .true.;  return
      endif
    else

      ! The slices are independent (disjoint data, one private bunch scratch per
      ! thread), so the loop runs slice-parallel; track1_bunch's own particle-level
      ! OMP region nests inside and, with nesting off (the OpenMP default), runs
      ! serial per thread -- coarse slice granularity replaces fine particle
      ! granularity, and each slice's arithmetic is untouched (bit-for-bit; the
      ! thread-identity checks cover it). Radiation FLUCTUATIONS draw from the one
      ! shared RNG stream inside track1, whose draw ORDER must stay fixed: that
      ! (rare, check-mode) configuration keeps the serial loop.

      if (bmad_com%radiation_fluctuations_on) then
        do is = 1, nslice
          call fel_slice_to_bunch (fbeam, fbeam%slice(is), ele, bunch, err)
          if (err) then
            err_flag = .true.;  return
          endif
          call track1_bunch (bunch, ele, err)
          if (err) then
            call out_io (s_error$, r_name, 'TRACKING ERROR IN ELEMENT: ' // trim(ele%name))
            err_flag = .true.;  return
          endif
          call fel_bunch_to_slice (bunch, ele, fbeam%slice(is), err)
          if (err) then
            err_flag = .true.;  return
          endif
        enddo
      else
        err = .false.
        !$OMP parallel do firstprivate(err) reduction(.or.: err_flag) schedule(static)
        do is = 1, nslice
          if (.not. err) then

            ! The scratch bunch is BLOCK-LOCAL: freshly default-initialized each
            ! iteration by Fortran's own semantics. (An OMP private clause on the
            ! subroutine-level scratch left the derived type's components
            ! improperly initialized under gfortran -- found as a deterministic
            ! 2x bunching shift -- so the scratch's definition lives here, where
            ! no clause semantics are involved.)

            block
              type (bunch_struct) bunch_l
              call fel_slice_to_bunch (fbeam, fbeam%slice(is), ele, bunch_l, err)
              if (.not. err) call track1_bunch (bunch_l, ele, err)
              if (.not. err) call fel_bunch_to_slice (bunch_l, ele, fbeam%slice(is), err)
              if (err) err_flag = .true.
            end block
          endif
        enddo
        !$OMP end parallel do
        if (err_flag) then
          call out_io (s_error$, r_name, 'TRACKING ERROR IN ELEMENT: ' // trim(ele%name))
          return
        endif
      endif
    endif

    fbeam%phi0 = fbeam%phi0 + ele%value(l$) * &
                    fel_phi0_rate(ks, ks * 0.5_rp / gamma0_ref**2, fel_p0_mc(fbeam))

    do ih = 1, n_harm      ! Each field diffracts at its own wavelength; through a
                         ! geometry break the light goes the CHORD, not the arc, and
                         ! the correction lands on the break's last element.
    call wavefront_drift (ffield(ih)%wf, ele%value(l$) - light_corr(ie), err)
    if (err) exit
  enddo

  ! Absolute-time phasing (manual sec:phasing; bmad_com's global switch through
  ! Bmad's own resolver): keep the real beam-vs-light carrier phase of this break --
  ! the drift slip plus any geometric (chicane) delay -- where the relative mode
  ! re-anchors. Whole turns wrap; no floors needed on a phase.

  if (absolute_time_tracking(ele)) then
    fbeam%phi0 = fbeam%phi0 - phase_rate * ele%value(l$) - twopi * light_corr(ie) / wf%wavelength
  endif
    if (err) then
      err_flag = .true.;  return
    endif
    ! The chamber does not end where the undulator does: the wake's energy loss applies
    ! through seam interludes too, as one kick of the element's length (Genesis applies
    ! it every step, and an interlude is one step).

    call fel_wake_apply (coll%wake, fbeam, ele%value(l$))

    z_now = z_now + ele%value(l$)
    call apply_slippage_banked (ele_slip(ie))
    if (err_flag) return

    call do_migrate ()
    if (err_flag) return
    if (fel_comb_take(comb, z_now, run%z_last_rec, .true.)) then
      call take_stats_record (.true.)
      if (err_flag) return
      call write_diag_rows()
    endif
    call progress_line (.true., 1, 1)
    call end_of_element ()
    if (err_flag) return

  else

    ! Genesis's own interlude model, transcribed, for pricing what the seam changes.

    qf = 0
    if (ele%key == quadrupole$) qf = ele%value(k1$)
    call fel_track_interlude_genesis (qf, ele%value(l$), fbeam, ffield, coll, err)
    if (absolute_time_tracking(ele)) then      ! Absolute-time phasing (sec:phasing).
      fbeam%phi0 = fbeam%phi0 - phase_rate * ele%value(l$)
    endif
    if (associated(wake_src)) then
      call apply_bmad_wake_kick (wake_src)
      if (err_flag) return
    endif
    if (err) then
      err_flag = .true.;  return
    endif
    z_now = z_now + ele%value(l$)
    call apply_slippage_banked (ele_slip(ie))
    if (err_flag) return

    call do_migrate ()
    if (err_flag) return
    if (fel_comb_take(comb, z_now, run%z_last_rec, .true.)) then
      call take_stats_record (.true.)
      if (err_flag) return
      call write_diag_rows()
    endif
    call progress_line (.true., 1, 1)
    call end_of_element ()
    if (err_flag) return
  endif
enddo

if (write_diag) close (iu_diag)
if (any_unavg) close (iu_ledger)
if (run%coll%wake%on) close (iu_wake)

if (migrate) then
  call out_io (s_info$, r_name, 'Migration moved \i0\ particles; dropped charge \es12.4\ C off the window ends.', &
               i_array = [n_moved_tot], r_array = [charge_dropped_tot])
  if (migrate_check) then
    call out_io (s_info$, r_name, 'Worst whole-beam bunching deviation across migrations: \es10.2\ ', &
                 r_array = [b_dev_max])
  endif
endif


!------------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!+
! Subroutine apply_bmad_wake_kick (wake_ele)
!
! Routine to apply one whole-window application of an element's Bmad sr wake, as a
! pure kick (manual sec:seamwake): concatenate the slices into global window
! coordinates, let Bmad's own machinery order and kick (track1_sr_wake: pseudomode
! accumulation head to tail, z_long binned FFT), split back holding theta --
! Genesis's convention for wake energy loss, the same z rescale fel_wake_apply_slice
! does -- so the phase every deposition sees is continuous through the kick. Used
! inside wigglers (mid-element) and after genesis-model interludes; Bmad-model
! interludes instead go through track1_bunch, where the wake applies at ds_wake in
! Bmad's own chart.
!-

subroutine apply_bmad_wake_kick (wake_ele)

type (ele_struct) wake_ele
logical err_w

!

call fel_concat_slices (fbeam, wake_ele, run%wake_bunch, run%wake_beta0, err_w)
if (err_w) then
  err_flag = .true.;  return
endif
call order_particles_in_z (run%wake_bunch)
call track1_sr_wake (run%wake_bunch, wake_ele)
call fel_split_slices (run%wake_bunch, wake_ele, fbeam, run%wake_beta0, .true., err_w)
if (err_w) then
  err_flag = .true.;  return
endif
end subroutine apply_bmad_wake_kick

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine do_migrate ()
!
! Routine to run slice migration at the per-element stride, serial, between the
! parallel regions (the thread check stays untouched). Called AFTER z_now is advanced,
! so per-event drop reports carry the z of the diagnostic record they precede -- the
! conservation timeline reconstructs exactly from the log. With migrate_check, the
! whole-beam weighted phasor S = sum(w e^{i theta}) must satisfy
! S_before = S_after + S_dropped to rounding: every mover's phase shifts by an exact
! multiple of 2*pi*sample and a drop removes exactly its own term, so any deviation
! beyond rounding is a bookkeeping bug (wrong z adjustment, weight not moved), not
! statistics.
!-

subroutine do_migrate ()

real(rp) chd, sb_re, sb_im, sa_re, sa_im, d_re, d_im, wsum
integer nm

if (.not. migrate) return

if (migrate_check) call whole_beam_phasor (sb_re, sb_im, wsum)

call fel_migrate_slices (fbeam, ks, nm, chd, d_re, d_im, err)
if (err) then
  err_flag = .true.;  return
endif
n_moved_tot = n_moved_tot + nm
charge_dropped_tot = charge_dropped_tot + chd

! Migration changes the current profile, which the wake convolution was hoisted on
! (the hoist predates migration): recompute at this stride. Every recompute appends
! a z-stamped block to <out_root>.wake.txt, so "the wake followed the currents"
! is a structural fact a check can parse without reimplementing the convolution.

if (nm > 0 .and. coll%wake%on) then
  call fel_wake_update (coll%wake, fbeam)
  call fel_write_wake_block (run, z_now)
endif
if (chd > 0) then
  call out_io (s_info$, r_name, 'Migration dropped \es22.14\ C off the window ends at z = \es22.14\ m.', &
               r_array = [chd, z_now])
endif

if (migrate_check .and. (nm > 0 .or. chd > 0)) then
  call whole_beam_phasor (sa_re, sa_im, wsum)
  if (wsum > 0) then
    b_dev_max = max(b_dev_max, sqrt((sb_re - sa_re - d_re)**2 + (sb_im - sa_im - d_im)**2) / wsum)
  endif
endif

end subroutine do_migrate

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine whole_beam_phasor (s_re, s_im, wsum)
!
! Routine to compute the whole-beam weighted phasor sum(w e^{i theta}) and the total
! weight over all slices. Exactly conserved across migration (moves shift phases by
! 2*pi*sample multiples; drops are accounted separately), which is what migrate_check
! verifies.
!-

subroutine whole_beam_phasor (s_re, s_im, wsum)

real(rp) s_re, s_im, wsum, theta, beta, p0_mc, w
integer is, ip

s_re = 0; s_im = 0; wsum = 0
p0_mc = fel_p0_mc(fbeam)

do is = 1, size(fbeam%slice)
  do ip = 1, fbeam%slice(is)%n
    w = fbeam%slice(is)%weight(ip)
    beta = fel_beta_of(p0_mc, fbeam%slice(is)%pz(ip))
    theta = fbeam%phi0 + ks * fbeam%slice(is)%z(ip) / beta
    s_re = s_re + w * cos(theta)
    s_im = s_im + w * sin(theta)
    wsum = wsum + w
  enddo
enddo

end subroutine whole_beam_phasor

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine write_diag_rows ()
!
! Routine to write one diag row per slice, slices in time-window order: beam slice is
! against field slice fel_field_index(slip, is, nslice), the unrotation of manual
! sec:slippage. The values are the ones the stats loop just evaluated with the SAME
! fel_field_diag and fel_slice_diag calls, slice-parallel (each slice's arithmetic
! identical to the old serial sweep, so this file is bit-for-bit what it always was);
! this routine only prints. take_stats_record must have run for this record first.
!-

subroutine write_diag_rows ()

integer is

if (.not. write_diag) return

do is = 1, nslice
  write (iu_diag, '(es24.16, i8, 10es24.16)') z_now, is, fpow_arr(is), fonax_arr(is), &
        bdiag_arr(is)%bunching, bdiag_arr(is)%bunching_phase, bdiag_arr(is)%mean_energy, &
        bdiag_arr(is)%sigma_energy, bdiag_arr(is)%sigma_x, bdiag_arr(is)%sigma_y, &
        bdiag_arr(is)%current, bdiag_arr(is)%n_eff
enddo

end subroutine write_diag_rows

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine write_ledger_row ()
!
! Routine to write one energy-ledger row (unaveraged mode): beam energy RELATIVE to
! the reference, sum(w*(gamma-gamma0))*me [C*eV = J] -- relative so the per-record
! change is not differenced off a large baseline at its own summation-rounding floor
! (manual sec:numerics) -- total window field energy sum(P_is)*slice_spacing/c [J],
! and the kick-side change dE_step returned by fel_unavg_step. The last columns are the
! cumulative energy transmitted out of the window by slippage (banked at the zero fill
! in fel_apply_slippage), the cumulative spontaneous deposit energy sum|dE_src|^2 (the
! one field-energy term the kick/deposit duality does not charge to the beam), and the
! cumulative energy the beam radiated away under bmad_com's radiation switches (the
! ACTUAL drawn sums, not expectations, so closure stays exact; zero with the switches
! off). In a time-dependent run the window is an open system, and the EXACTLY closing
! quantity is E_beam + U_field + U_escaped - U_spont + E_radiated. Wakes would be a
! second, unbanked exit channel; this ledger only exists where they are refused.
!-

subroutine write_ledger_row ()

real(rp) e_beam, u_field, p_mc_l, g0_l
integer is, ip

! The field power per slice comes from the stats loop's fel_field_diag evaluation
! (same call, same unrotation, same summation order over slices -- bit-identical);
! take_stats_record must have run for this record first.

e_beam = 0
u_field = 0
g0_l = fel_gamma0(fbeam)
do is = 1, nslice
  do ip = 1, fbeam%slice(is)%n
    p_mc_l = fel_p0_mc(fbeam) * (1 + fbeam%slice(is)%pz(ip))
    e_beam = e_beam + fbeam%slice(is)%weight(ip) * (sqrt(p_mc_l**2 + 1) - g0_l) * m_electron
  enddo
  u_field = u_field + fpow_arr(is) * fbeam%slice_spacing / c_light
enddo

write (iu_ledger, '(7es24.16)') z_now, e_beam, u_field, dE_step, slip%u_escaped, u_spont_cum, e_rad_cum

end subroutine write_ledger_row

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine apply_radiation ()
!
! Routine to apply spontaneous radiation inside FEL elements, honoring Bmad's GLOBAL
! switches bmad_com%radiation_damping_on / %radiation_fluctuations_on -- the same
! switches every Bmad tracking path honors (interludes get theirs through track1; this
! covers the custom-tracked FEL step, whose radiation reaction cannot emerge from its
! Newton-Lorentz equations of motion). Per record:
!
!   damping:      dgamma_j = -(2/3) r_e gamma_j^2 ku^2 aw^2 * INT g^2 ds
!                 (each particle's own gamma; the envelope integral over the record, so
!                 the unaveraged ramps radiate by their actual strength; averaged g = 1)
!   fluctuations: the standard undulator quantum-diffusion form -- variance
!                 1.015e-27 * ku^3 aw^2 F(aw) gamma0^4 per meter with Genesis's
!                 Incoherent.cpp F(aw) fits (the Saldin closed form) -- drawn GAUSSIAN.
!                 Genesis draws uniform scaled by sqrt(3) to reach this SAME variance;
!                 the sqrt(3) normalizes the uniform draw and must not appear with a
!                 Gaussian (the physical limit, a merit choice). ONE DRAW PER BEAMLET,
!                 exactly as
!                 Genesis: independent per-particle kicks would break the quiet start's
!                 per-beamlet harmonic cancellation. Draws happen SERIALLY in fixed
!                 slice order from the one seeded stream, so results are independent of
!                 thread count; only the application parallelizes.
!
! The actual drawn energy accumulates into e_rad_cum (the ledger's E_radiated column):
! drawn sums, not expectations, so the TD closure stays exact. In the unaveraged mode
! this double-counts the grid-captured band (~3% at the reference grid, bounded live by
! check_spontaneous.py) -- accepted, smaller than the angular-distribution estimate's
! own uncertainty.
!-

subroutine apply_radiation ()

real(rp) c_loss, fform, sig, gam, p_mc_r, dgam, intg2, s0, g_env, gp_env
integer isl, ipr, ibl, nb, nbl_max
!

if (.not. (bmad_com%radiation_damping_on .or. bmad_com%radiation_fluctuations_on)) return

! The envelope integral over this record: the substep-grid midpoint sum for the
! unaveraged mode (matching its own integration grid), dz exactly for the averaged.

if (fel_mode(ie) == fel_unaveraged$) then
  intg2 = 0
  do ipr = 1, run%ustate%nsub
    s0 = (run%ustate%s - und%dz) + (ipr - 0.5_rp) * run%ustate%dsub
    g_env = fel_unavg_envelope(run%ustate, s0, gp_env)
    intg2 = intg2 + g_env**2 * run%ustate%dsub
  enddo
else
  intg2 = und%dz
endif

c_loss = (2.0_rp / 3.0_rp) * r_e * und%ku**2 * und%aw**2
if (und%helical) then
  fform = 1.42_rp * und%aw + 1 / (1 + 1.5_rp * und%aw + 0.95_rp * und%aw**2)
else
  fform = 1.697_rp * und%aw + 1 / (1 + 1.88_rp * und%aw + 0.8_rp * und%aw**2)
endif
sig = 0
if (bmad_com%radiation_fluctuations_on) then
  sig = sqrt(1.015e-27_rp * und%ku**3 * und%aw**2 * fform * gamma0_ref**4 * intg2)
endif

nb = max(1, fbeam%nbins)
if (.not. allocated(e_rad_slice)) allocate (e_rad_slice(nslice))
e_rad_slice = 0

if (bmad_com%radiation_fluctuations_on) then
  nbl_max = 0
  do isl = 1, nslice
    nbl_max = max(nbl_max, (fbeam%slice(isl)%n + nb - 1) / nb)
  enddo
  if (.not. allocated(rad_kick)) allocate (rad_kick(nbl_max, nslice))
  do isl = 1, nslice
    do ibl = 1, (fbeam%slice(isl)%n + nb - 1) / nb
      call ran_gauss (rad_kick(ibl, isl))
    enddo
  enddo
endif

!$OMP parallel do private(ipr, p_mc_r, gam, dgam)
do isl = 1, nslice
  do ipr = 1, fbeam%slice(isl)%n
    p_mc_r = fel_p0_mc(fbeam) * (1 + fbeam%slice(isl)%pz(ipr))
    gam = sqrt(p_mc_r**2 + 1)
    dgam = 0
    if (bmad_com%radiation_damping_on) dgam = dgam - c_loss * gam**2 * intg2
    if (bmad_com%radiation_fluctuations_on) dgam = dgam + sig * rad_kick((ipr-1)/nb + 1, isl)
    e_rad_slice(isl) = e_rad_slice(isl) - fbeam%slice(isl)%weight(ipr) * dgam * m_electron
    gam = gam + dgam
    fbeam%slice(isl)%pz(ipr) = (sqrt(gam**2 - 1) - fel_p0_mc(fbeam)) / fel_p0_mc(fbeam)
  enddo
enddo
!$OMP end parallel do

e_rad_cum = e_rad_cum + sum(e_rad_slice)

end subroutine apply_radiation

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine progress_line (at_element_end, i_step, n_step)
!
! Routine to print one progress line to stdout: where the walk is and what the light
! and beam are doing, so the slow modes (the unaveraged mode runs ~30x the averaged)
! show signs of life. Element boundaries always print; inside elements a wall-clock
! throttle (2 s) keeps fast runs quiet. All numbers are read from the stats row just
! taken.
!-

subroutine progress_line (at_element_end, i_step, n_step)

logical at_element_end
integer i_step, n_step
integer(8) now
real(rp) elapsed
character(200) line

call system_clock (now)
if (.not. at_element_end .and. real(now - prog_count_last, rp) / prog_rate < 2.0_rp) return
prog_count_last = now
elapsed = real(now - prog_count0, rp) / prog_rate

! With no stats record to read (comb_ds_save < 0), the walk still shows signs of
! life: the same line without the physics numbers the records would carry.

if (stats%irec == 0) then
  write (line, '(a, f5.1, a, f8.3, a, i0, a, i0, 3a, i0, a, i0, a, i0, a)') &
        'progress: ', 100 * z_now / lat_length, '%  z = ', z_now, ' m  ele ', ie, '/', &
        branch%n_ele_track, ' ', trim(ele%name), '  step ', i_step, '/', n_step, &
        '  t = ', nint(elapsed), ' s'
else
  write (line, '(a, f5.1, a, f8.3, a, i0, a, i0, 3a, i0, a, i0, a, es9.2, a, es9.2, a, f8.5, a, i0, a)') &
        'progress: ', 100 * z_now / lat_length, '%  z = ', z_now, ' m  ele ', ie, '/', &
        branch%n_ele_track, ' ', trim(ele%name), '  step ', i_step, '/', n_step, &
        '  P = ', sum(stats%f_power(:, stats%irec)), ' W  U = ', sum(stats%f_energy(:, stats%irec)), &
        ' J  <|b|> = ', sum(stats%bunching(:, stats%irec)) / nslice, '  t = ', nint(elapsed), ' s'
endif
call out_io (s_blank$, r_name, trim(line))

end subroutine progress_line

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine take_stats_record (with_angles)
!
! Routine to take one stats record at the current z_now through fel_stats_record.
! with_angles fills the field theta moments (element ends). Sets the host err_flag
! on error.
!-

subroutine take_stats_record (with_angles)

logical with_angles, serr

call fel_stats_record (stats, fbeam, ffield, z_now, with_angles, bdiag_arr, fpow_arr, fonax_arr, serr)
if (serr) then
  err_flag = .true.;  return
endif
end subroutine take_stats_record

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine end_of_element ()
!
! Routine to handle the end of an element: the evaluated bunch_params (the Tao
! end-of-element pattern) and any requested dumps. Mid-run field dumps unrotate
! exactly as the final dump does -- the rotation is a gauge, and re-zeroing slip%first
! keeps the mapping consistent.
!-

subroutine end_of_element ()

logical eerr
character(500) fname

!

call fel_stats_element_end (stats, fbeam, ffield, ele, z_now, eerr)
if (eerr) then
  err_flag = .true.;  return
endif
if (dump_beam_here(ie)) then
  write (fname, '(2a, i0, 3a)') trim(out_root), '-at', ie, '-', trim(ele%name), '.par.h5'
  call fel_write_genesis4_beam (fbeam, trim(fname), eerr)
  if (eerr) then
    err_flag = .true.;  return
  endif
endif

if (dump_field_here(ie)) then
  write (fname, '(2a, i0, 2a)') trim(out_root), '-at', ie, '-', trim(ele%name)
  call fel_dump_field_set (run, trim(fname), eerr)
  if (eerr) then
    err_flag = .true.;  return
  endif
endif

end subroutine end_of_element

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine apply_slippage_banked (slippage)
!
! Routine to apply slippage to every field of the set and drain the escape bank when
! the escaped field is kept. Every field slips by the same amount in
! fundamental-wavelength units (one window, lockstep rotation; the harm argument only
! fixes the escape bank's slice light-time for a harmonic field).
!-

subroutine apply_slippage_banked (slippage)

real(rp) slippage
integer ihh
logical err_b

!

if (keep_escaped_field) then
  do ihh = 1, n_harm
    call fel_apply_slippage (ffield(ihh)%slip, ffield(ihh)%wf, slippage, ffield(ihh)%bank, &
                             ffield(ihh)%harm)
    call fel_drain_bank (run, ihh, err_b)
    if (err_b) then
      err_flag = .true.;  return
    endif
  enddo
else
  do ihh = 1, n_harm
    call fel_apply_slippage (ffield(ihh)%slip, ffield(ihh)%wf, slippage, harm = ffield(ihh)%harm)
  enddo
endif

end subroutine apply_slippage_banked

end subroutine track_fel_line

end module fel_track_line_mod
