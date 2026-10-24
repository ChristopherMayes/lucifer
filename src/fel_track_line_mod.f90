!+
! Module fel_track_line_mod
!
! The walk. track_fel_line(run) is the main loop of the FEL tracker.
! FEL segments are stepped with the transcribed Genesis physics. Every other element
! goes through the Bmad seam, or through the transcribed Genesis interlude model.
! Slippage and phasing follow the precomputed schedule. Stats and diag records are
! taken at Genesis's record positions.
!
! The routine can be called repeatedly in one process. All state lives in run and the
! RNG is seeded at init. Errors return through err_flag. Nothing here stops.
!-

module fel_track_line_mod

use fel_struct
use fel_io_mod
use fel_timer_mod
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
! Routine to walk the lattice from run%i_start through run%i_end, the resolved
! tracking window.
! The schedule was built on the full lattice, so a windowed run composes exactly with
! the full one. The routine opens and closes the run's diag, ledger and wake streams.
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
integer iu_mig

type (bunch_struct) bunch
type (fel_und_struct) und
real(rp), allocatable :: e_rad_slice(:), rad_kick(:,:)
real(rp) gamma0_ref, phase_rate, ks, qf, und_slip_step, lat_length, dE_step, dU_step, comb
integer nslice, n_harm, ie, is, ih, istep
integer(8) prog_count0, prog_count_last, prog_rate
real(rp) prog_power, prog_energy, prog_bunch
logical prog_header_done
logical write_diag, keep_escaped_field, migrate, migrate_check, any_unavg, two_pol, err
character(400) out_root
character(16) interlude_model
type (ele_struct), target :: slice_ele          ! One interlude piece, element_slice_iterator's.
type (ele_struct), pointer :: ele_p             ! The element the piece's arithmetic reads.
integer i_piece, n_piece
real(rp) l_piece, f_piece, z_ele0
character(*), parameter :: r_name = 'track_fel_line'

!

err_flag = .false.
call fel_tic (fel_t_walk$)
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

! The unaveraged energy ledger (fel-physics.md sec-unaveraged). One row per record
! step inside FEL segments: total weighted beam energy, total window field energy,
! and the kick-side energy change of the step. The ledger check holds
! d(E_beam + U_field) to the floor in docs/validation.md.

if (any_unavg) then
  open (newunit = iu_ledger, file = trim(out_root) // '.ledger.txt', action = 'write')
  write (iu_ledger, '(a)') '#         z          E_beam_rel [J]          U_field [J]           dE_kick [J]          U_escaped [J]          U_spont [J]         E_radiated [J]'
endif

n_moved_tot = 0
charge_dropped_tot = 0
b_dev_max = 0

! The migration stream, one row per event, opened like the diag and ledger streams.

if (migrate) then
  open (newunit = iu_mig, file = trim(out_root) // '.migration.txt', action = 'write')
  write (iu_mig, '(a)') '# Slice migration, one row per event. Machine-readable; stdout is not.'
  write (iu_mig, '(a)') '#          s [m]        moved      charge_dropped [C]     phasor_deviation' // &
                        '        min_margin [rad]'
endif

! Progress goes to stdout, throttled by wall clock. Slow modes print a steady trickle
! and fast runs print only the element boundaries. The numbers come from the stats row
! just taken, so nothing extra is computed.

call system_clock (prog_count0, prog_rate)
prog_count_last = prog_count0
prog_power = 0;  prog_energy = 0;  prog_bunch = 0
prog_header_done = .false.
lat_length = branch%ele(branch%n_ele_track)%s

z_now = branch%ele(run%i_start - 1)%s      ! 0 for a full run, else the window's entry face.
call set_field_position ()
ie = run%i_start - 1                      ! The entry face, element 0 for a full run. The
                                          !   initial record needs ie before the loop sets it.
if (fel_comb_take(comb, z_now, run%z_last_rec, .false.)) then
  call take_stats_record (.true.)   ! Evaluates the diag rows too. The writer prints.
  if (err_flag) return
  call write_diag_rows()            ! Initial record, matching Genesis's diag before the first step.
endif



do ie = run%i_start, run%i_end
  ele => branch%ele(ie)

  ! Zero length elements (Bmad's end marker, for one) get no step and no diagnostic
  ! record. Genesis's unrolled lattice has no counterpart for them. The exception is a
  ! wake applying there (a zero-length wake element is a standard Bmad idiom). Such an
  ! element takes the interlude wake path below.
  !
  ! wake_src resolves the wake through lords. A wake on a superimposed or split element
  ! lives on the lord, and ele%wake is null on its slaves (pointer_to_wake_ele). The
  ! resolution picks exactly one slave of a split lord, the one containing the lord's
  ! midpoint. So a split wake applies once, as in Bmad. Checking ele%wake directly would
  ! send lord wakes to the per-slice path, where Bmad applies them within single slices.
  !
  ! With bmad_com%sr_wakes_on off, track1_bunch applies no wake. The whole-window
  ! concatenation below then has only its serial cost, so a wake-carrying element routes
  ! like any other interlude, through the parallel slice loop. On a line of many wake
  ! pipes the switch decides between the serial and parallel interlude paths for the
  ! whole lattice (docs/performance.md).
  !
  ! A zero-length element that cannot act is skipped. That is a marker, or a pipe whose
  ! wake scales with its length and so kicks by zero. A wake assigned over a range or a
  ! class lands on such pipes routinely. Tracking one would cost the serial interlude
  ! path, a record, an element end and a repeated coords/s, for no change. A zero-length
  ! element that can act (a thin kicker, a multipole, a long-range wake) is tracked.
  ! Both decisions are fel_struct's, and the setup's record-count precompute makes the
  ! same two from the same functions.

  wake_src => fel_interlude_wake (ele)
  if (fel_element_inert (ele, wake_src)) cycle

  ! Space charge is this element's own attribute, resolved at setup into run%sc_here
  ! with Bmad's master switch folded in. The solvers read coll%efield%active and zero
  ! their outputs when it is false. So an element without space charge tracks the
  ! space-charge-free path bit for bit.

  coll%efield%active = run%sc_here(ie)

  if (is_fel(ie)) then

    ! FEL segment: Genesis's unroll in the element's own step. The step comes from
    ! Bmad's ds_step and num_steps attributes, whose bookkeeper computes Genesis's
    ! num_steps = round(l/ds_step) (attribute_bookkeeper.f90). There is no namelist step
    ! size. The steps are equal, and slippage follows each step's field solve. Any
    ! end-of-lattice fixup lands on the last step.

    und = und_of(ie)
    ! Transverse transport in an averaged element uses Bmad's own kernel by default.
    ! transport_model = 'genesis' selects Genesis4's maps, for comparison runs. The
    ! unaveraged mode integrates the field and has no map to choose.
    und%bmad_transport = (run%global%transport_model == 'bmad')
    und_slip_step = (1 + und%aw**2) / (2 * gamma0_ref**2 * wf%wavelength)
    und%nstep = max(1, nint(ele%value(num_steps$)))
    und%dz = ele%value(l$) / und%nstep

    if (fel_mode(ie) == unaveraged$) then
      ! The concatenated wake kick would meet the quiver-carrying chart mid-segment.
      ! Nothing in this mode needs element wakes, so the run is refused.
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

    ! The off-phase knob (fel-physics.md sec-phasing). A displaced element sees the extra
    ! upstream-break phase at entry. It gives the phase back at its last step, before the
    ! migration there, since the downstream break is shorter by the same delta. So the
    ! anchor stays nominal downstream. Sign: positive z_offset = a longer upstream break
    ! = more beam delay = theta backwards. This is Genesis's phase-shifter convention,
    ! anchored by the cross-code phi scan.

    if (fel_zoff(ie) /= 0) fbeam%phi0 = fbeam%phi0 - phase_rate * fel_zoff(ie)

    ! Device residency (fel_device_mod). With the backend on and fp32_check off, the
    ! beam and the field upload here and stay resident through the element. Element
    ! wakes are refused here, since only this walk resolves them through the lords.

    if (run%dev%on .and. .not. run%fp32%on .and. fel_mode(ie) /= unaveraged$) then
      if (associated(wake_src)) then
        call out_io (s_error$, r_name, 'DEVICE = "' // trim(run%global%device) // &
                     '" DOES NOT COVER ELEMENT WAKES. AT ELEMENT: ' // trim(ele%name))
        err_flag = .true.;  return
      endif
      call fel_device_element_begin (run%dev, fbeam, ffield)
    endif

    do istep = 1, und%nstep
      if (fel_mode(ie) == unaveraged$) then
        call fel_tic (fel_t_unavg$)
        call fel_unavg_step (und, run%ustate, fbeam, run%ffield(1:1), und%dz, istep == 1, &
                             istep == und%nstep, dE_step, dU_step, run%fp32, run%dev, err)
        call fel_toc (fel_t_unavg$)
        u_spont_cum = u_spont_cum + dU_step
      else
        call fel_track_und_step (und, fbeam, ffield, coll, err, run%fp32, run%dev)
      endif
      if (err) then
        err_flag = .true.;  return
      endif
      call apply_radiation ()

      ! Element sr wake, once per passage as in Bmad. One kick at the step nearest
      ! mid-element, scaled to the full element length (scale_with_length uses ele's l),
      ! across the whole window (fel-physics.md sec-seamwake). It is a direct kick with
      ! no transport, since this walk owns transport inside wigglers.

      if (associated(wake_src) .and. istep == (und%nstep + 1)/2) then
        call apply_bmad_wake_kick (wake_src)
        if (err_flag) return
      endif

      z_now = z_now + und%dz
      call set_field_position ()
      if (istep == und%nstep) then
        call apply_slippage_banked (und%dz * und_slip_step + ele_slip(ie))
        if (err_flag) return
      else
        call apply_slippage_banked (und%dz * und_slip_step)
        if (err_flag) return
      endif
      ! The offset is given back before the migration re-slices the beam. Migration
      ! reads each particle's phase against its slice's window. With the displaced phase
      ! still applied, every mover lands outside its new slice once the nominal phase
      ! returns (docs/validation.md).

      if (istep == und%nstep) then
        if (fel_zoff(ie) /= 0) fbeam%phi0 = fbeam%phi0 + phase_rate * fel_zoff(ie)
        call do_migrate ()
        if (err_flag) return
      endif
      if (fel_comb_take(comb, z_now, run%z_last_rec, istep == und%nstep)) then

        ! At the comb's positions the resident state comes back to the host. The stats,
        ! the diag rows and everything downstream read the refreshed FP64 arrays.
        ! Between rows the host state is stale. Off the device this does nothing.

        call fel_device_readback (run%dev, fbeam, ffield)
        if (run%dev%dep_breach) then;  err_flag = .true.;  return;  endif
        call take_stats_record (istep == und%nstep)
        if (err_flag) return
        if (fel_mode(ie) == unaveraged$) call write_ledger_row ()
        call write_diag_rows()
      endif
      if (istep /= und%nstep) call progress_line (.false., istep, und%nstep)
    enddo
    call fel_device_element_end (run%dev, fbeam, ffield)   ! The residency boundary: z re-forms in FP64.

    ! A deposit met a particle below the gamma its bound assumes and dropped it. So this
    ! element's field is wrong, and nothing downstream may read it. The run stops here,
    ! before a stats row, a dump or an exit line is written from it.

    if (run%dev%dep_breach) then;  err_flag = .true.;  return;  endif
    call end_of_element ()              ! Fills the element-end row before this element's
    if (err_flag) return                !   last progress row reads it.
    call progress_line (.true., und%nstep, und%nstep)

  elseif (interlude_model == 'bmad') then

    ! The interlude's own steps. Tao's comb has two halves. comb_ds_save carries the
    ! first: save_a_bunch_step selects among positions a tracker already reaches. This is
    ! the second, tao_lattice_calc_mod's n_slice branch. It cuts an element into pieces
    ! with element_slice_iterator where the tracker would cross it in one map. Without
    ! it a drift reaches one position, its end, and no comb setting adds another.
    !
    ! Tao does not cut where CSR or space charge would see the cut. The same holds here
    ! for the element short-range wake. That wake is Bmad's once-per-passage kick of the
    ! element's length and would become one kick a piece. The chamber wake is an energy
    ! loss proportional to the length, applied additively in gamma. So its pieces sum to
    ! the element's own kick, and it does not stop the cut.

    n_piece = 1
    if (.not. associated(wake_src)) &
            n_piece = fel_interlude_pieces(run%global%interlude_ds_step, ele%value(l$))
    z_ele0 = z_now

    do i_piece = 1, n_piece
      if (n_piece == 1) then
        ele_p => ele
      else
        call element_slice_iterator (ele, branch%param, i_piece, n_piece, slice_ele)
        ele_p => slice_ele
      endif
      l_piece = ele_p%value(l$)

      ! The break's chord-vs-arc correction belongs to the whole break and lands on its
      ! last element (fel_setup_mod's close_geometry_break). Here it lands on that
      ! element's last piece, so the pieces sum to the element. The slippage is
      ! proportional to length and splits.

      f_piece = 0
      if (i_piece == n_piece) f_piece = 1

      ! The seam. Bmad tracks each slice's bunch (coordinate copies in and out).
      ! wavefront_drift moves the field (every slice, rotation-invariant). The common
      ! phase phi0 advances by the reference rate, with Genesis's drift surrogate
      ! ks/(2*gamma0^2) as the reference wavenumber.

      call fel_tic (fel_t_seam$)
      if (associated(wake_src)) then

        ! Wake-carrying interlude: all slices as one bunch in global window coordinates,
        ! through Bmad's track1_bunch. track1_bunch applies the sr wake at ds_wake with
        ! the whole window visible head to tail (fel-physics.md sec-seamwake). Every
        ! other element takes the per-slice path below.

        call fel_concat_slices (fbeam, ele_p, run%wake_bunch, run%wake_beta0, err)
        if (err) then
          err_flag = .true.;  return
        endif
        call track1_bunch (run%wake_bunch, ele_p, err)
        if (err) then
          call out_io (s_error$, r_name, 'TRACKING ERROR IN ELEMENT: ' // trim(ele%name))
          err_flag = .true.;  return
        endif
        call fel_split_slices (run%wake_bunch, ele_p, fbeam, run%wake_beta0, .false., err)
        if (err) then
          err_flag = .true.;  return
        endif
      else

        ! The slices are independent (disjoint data, one bunch scratch per thread), so
        ! the loop runs slice-parallel. track1_bunch's own particle-level OMP region
        ! nests inside. With nesting off (the OpenMP default) it runs serial per thread.
        ! Each slice's arithmetic is unchanged bit for bit, and the thread-identity
        ! checks cover it.
        !
        ! Radiation fluctuations draw from the one shared RNG stream inside track1.
        ! The draw order must stay fixed, so that configuration keeps the serial loop.

        if (bmad_com%radiation_fluctuations_on) then
          do is = 1, nslice
            call fel_slice_to_bunch (fbeam, fbeam%slice(is), ele_p, bunch, err)
            if (err) then
              err_flag = .true.;  return
            endif
            call track1_bunch (bunch, ele_p, err)
            if (err) then
              call out_io (s_error$, r_name, 'TRACKING ERROR IN ELEMENT: ' // trim(ele%name))
              err_flag = .true.;  return
            endif
            call fel_bunch_to_slice (bunch, ele_p, fbeam%slice(is), err)
            if (err) then
              err_flag = .true.;  return
            endif
          enddo
        else
          err = .false.
          !$OMP parallel do firstprivate(err) reduction(.or.: err_flag) schedule(static)
          do is = 1, nslice
            if (.not. err) then

              ! The scratch bunch is block-local, so Fortran default-initializes it
              ! each iteration.
              ! Note: An OMP private clause on a subroutine-level scratch leaves the
              ! derived type's components improperly initialized under gfortran. That
              ! gives a 2x bunching shift (docs/validation.md).

              block
                type (bunch_struct) bunch_l
                call fel_slice_to_bunch (fbeam, fbeam%slice(is), ele_p, bunch_l, err)
                if (.not. err) call track1_bunch (bunch_l, ele_p, err)
                if (.not. err) call fel_bunch_to_slice (bunch_l, ele_p, fbeam%slice(is), err)
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
      call fel_toc (fel_t_seam$)

      fbeam%phi0 = fbeam%phi0 + l_piece * &
                      fel_phi0_rate(ks, ks * 0.5_rp / gamma0_ref**2, fel_p0_mc(fbeam))

      call fel_tic (fel_t_drift$)
      do ih = 1, n_harm      ! Each field diffracts at its own wavelength. Through a
                             ! geometry break the light takes the chord and not the arc.
                             ! The correction lands on the break's last element.
        call wavefront_drift (ffield(ih)%wf, l_piece - light_corr(ie) * f_piece, err)
        if (err) exit
      enddo
      call fel_toc (fel_t_drift$)

      ! Absolute-time phasing (fel-physics.md sec-phasing), set by bmad_com's global switch
      ! through Bmad's own resolver. Keep the real beam-vs-light carrier phase of this
      ! break: the drift slip plus any geometric (chicane) delay. The relative mode
      ! re-anchors instead. Whole turns wrap, so a phase needs no floor.

      if (absolute_time_tracking(ele)) then
        fbeam%phi0 = fbeam%phi0 - phase_rate * l_piece - &
                     twopi * light_corr(ie) * f_piece / wf%wavelength
      endif
      if (err) then
        err_flag = .true.;  return
      endif
      ! The chamber does not end where the undulator does. The wake's energy loss applies
      ! through seam interludes too, as one kick of the element's length. Genesis applies
      ! it every step, and an interlude is one step.

      call fel_tic (fel_t_wake$)
      call fel_wake_apply (coll%wake, fbeam, l_piece)
      call fel_toc (fel_t_wake$)

      ! Several pieces take the piece boundary from the element's start. So the walk
      ! and the setup's record-count precompute reach the same z bit for bit, and the
      ! stats arrays stay exact-sized. An element's last boundary can then differ from
      ! z + l in its last bit. So one piece keeps the plain z + l form.

      if (n_piece == 1) then
        z_now = z_now + l_piece
      else
        z_now = z_ele0 + i_piece * (ele%value(l$) / n_piece)
      endif
      call set_field_position ()
      call apply_slippage_banked (ele_slip(ie) / n_piece)
      if (err_flag) return

      if (i_piece == n_piece) then
        call do_migrate ()
        if (err_flag) return
      endif
      if (fel_comb_take(comb, z_now, run%z_last_rec, i_piece == n_piece)) then
        call take_stats_record (i_piece == n_piece)
        if (err_flag) return
        call write_diag_rows()
      endif
      if (i_piece == n_piece) then
        call end_of_element ()            ! Fills the element-end row before the progress
        if (err_flag) return              !   row that reads it when there are no records.
      endif
      call progress_line (i_piece == n_piece, i_piece, n_piece)
      if (err_flag) return
    enddo

    ! element_slice_iterator leaves pointer components on its scratch element, and
    ! Fortran does not release those for a local variable. Each first call releases the
    ! previous element's. This releases the last element's.

    if (n_piece > 1) call deallocate_ele_pointers (slice_ele)

  else

    ! Genesis's own interlude model, transcribed, to compare against the seam.

    qf = 0
    if (ele%key == quadrupole$) qf = ele%value(k1$)
    call fel_tic (fel_t_gen_int$)
    call fel_track_interlude_genesis (qf, ele%value(l$), fbeam, ffield, coll, err)
    call fel_toc (fel_t_gen_int$)
    if (absolute_time_tracking(ele)) then      ! Absolute-time phasing (sec-phasing).
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
    call set_field_position ()
    call apply_slippage_banked (ele_slip(ie))
    if (err_flag) return

    call do_migrate ()
    if (err_flag) return
    if (fel_comb_take(comb, z_now, run%z_last_rec, .true.)) then
      call take_stats_record (.true.)
      if (err_flag) return
      call write_diag_rows()
    endif
    call end_of_element ()              ! Fills the element-end row before the progress
    if (err_flag) return                !   row that reads it when there are no records.
    call progress_line (.true., 1, 1)
    if (err_flag) return
  endif
enddo

if (write_diag) close (iu_diag)
if (any_unavg) close (iu_ledger)
if (run%coll%wake%on) close (iu_wake)
call fel_device_close_run (run%dev, run%fp32)   ! Its footer lines precede the stream's own.
call fel_fp32_close (run%fp32)

if (migrate) then
  write (iu_mig, '(a)') '#'
  write (iu_mig, '(a, i0)') 'moved ', n_moved_tot
  write (iu_mig, '(a, es24.15e3)') 'charge_dropped_total ', charge_dropped_tot
  if (migrate_check) write (iu_mig, '(a, es24.15e3)') 'worst_bunching_deviation ', b_dev_max
  close (iu_mig)
  call out_io (s_info$, r_name, 'Migration moved \i0\ particles, dropped \es12.4\ C off the window ends;' // &
               ' per-event rows in ' // trim(out_root) // '.migration.txt', &
               i_array = [n_moved_tot], r_array = [charge_dropped_tot])
endif

call fel_toc (fel_t_walk$)


!------------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!+
! Subroutine apply_bmad_wake_kick (wake_ele)
!
! Routine to apply an element's Bmad sr wake once to the whole window, as a pure kick
! (fel-physics.md sec-seamwake).
! The slices are concatenated into global window coordinates. Bmad's track1_sr_wake
! orders and kicks them (pseudomode accumulation head to tail, z_long binned FFT). The
! split back holds theta, Genesis's convention for wake energy loss and the same z
! rescale fel_wake_apply_slice does. So the phase every deposition sees is continuous
! through the kick.
! Used inside wigglers (mid-element) and after genesis-model interludes. Bmad-model
! interludes go through track1_bunch instead, where the wake applies at ds_wake in
! Bmad's own chart.
!-

subroutine apply_bmad_wake_kick (wake_ele)

type (ele_struct) wake_ele
logical err_w

!

call fel_tic (fel_t_wake$)
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
call fel_toc (fel_t_wake$)
end subroutine apply_bmad_wake_kick

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine do_migrate ()
!
! Routine to run slice migration at the per-element stride.
! Migration runs serially, between the parallel regions. It is called after z_now is
! advanced. So each drop report carries the z of the diagnostic record it precedes, and
! the conservation timeline can be rebuilt from the log.
!
! With migrate_check, the whole-beam weighted phasor S = sum(w e^{i theta}) must satisfy
!   S_before = S_after + S_dropped
! to rounding. Every mover's phase shifts by an exact multiple of 2*pi*sample, and a
! drop removes exactly its own term. Any larger deviation is a bookkeeping bug (wrong z
! adjustment, weight not moved). Statistics cannot cause it.
!-

subroutine do_migrate ()

real(rp) chd, sb_re, sb_im, sa_re, sa_im, d_re, d_im, wsum, dev, mrg
integer nm

if (.not. migrate) return
call fel_tic (fel_t_migration$)

! With the beam resident this is an element's last step. The migration re-slices the
! host arrays, so they must be current first, and no readback may follow. So the whole
! state comes back here, the field with it. Residency ends without the second readback
! fel_device_element_end would make. After this point the host holds the state, and
! the next element re-uploads. At interlude ends nothing is resident and this does
! nothing.

if (run%dev%resident) then
  call fel_device_readback (run%dev, fbeam, ffield)
  call fel_device_release (run%dev)
  if (run%dev%dep_breach) then;  err_flag = .true.;  return;  endif
endif

if (migrate_check) call whole_beam_phasor (sb_re, sb_im, wsum)

call fel_migrate_slices (fbeam, ks, nm, chd, d_re, d_im, mrg, err)
if (err) then
  err_flag = .true.;  return
endif
n_moved_tot = n_moved_tot + nm
charge_dropped_tot = charge_dropped_tot + chd

! Migration changes the current profile, and the wake convolution is computed from that
! profile. So the wake is recomputed at this stride. A particle dropped off the window's
! end changes the profile as a moved particle does. Every recompute appends a z-stamped
! block to <out_root>.wake.txt. A check can then confirm the wake followed the currents
! without redoing the convolution.

if ((nm > 0 .or. chd > 0) .and. coll%wake%on) then
  call fel_wake_update (coll%wake, fbeam)
  call fel_write_wake_block (run, z_now)
endif
dev = 0
if (migrate_check .and. (nm > 0 .or. chd > 0)) then
  call whole_beam_phasor (sa_re, sa_im, wsum)
  if (wsum > 0) then
    dev = sqrt((sb_re - sa_re - d_re)**2 + (sb_im - sa_im - d_im)**2) / wsum
    b_dev_max = max(b_dev_max, dev)
  endif
endif

! One row per event, at full precision, to the migration file (docs/user-guide.md): the
! charge that left the window ends and where, and the phase-continuity residual.
!
! The margin is the smallest phase distance any examined particle had to a slice
! boundary at this decision, dropped particles included. No frame written after the
! event can show it. A restart comparison reads it (docs/validation.md).

if (nm > 0 .or. chd > 0) then
  write (iu_mig, '(es22.14, i12, 3es24.15e3)') z_now, nm, chd, dev, mrg
endif
call fel_toc (fel_t_migration$)

end subroutine do_migrate

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine whole_beam_phasor (s_re, s_im, wsum)
!
! Routine to compute the whole-beam weighted phasor sum(w e^{i theta}) and the total
! weight over all slices.
! The phasor is conserved across migration. Moves shift phases by multiples of
! 2*pi*sample, and drops are counted separately. migrate_check verifies this.
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
! Routine to write one diag row per slice, slices in time-window order.
! Beam slice number is pairs with field slice fel_field_index(slip, is, nslice), the
! unrotation of fel-physics.md sec-slippage. The values are those the stats loop just
! evaluated, slice-parallel, with fel_field_diag and fel_slice_diag.
!
! Note: This routine only prints. take_stats_record must have run for this record first.
!-

subroutine write_diag_rows ()

integer is

if (.not. write_diag) return
call fel_tic (fel_t_stats$)

do is = 1, nslice
  write (iu_diag, '(es24.16, i8, 10es24.16)') z_now, is, fpow_arr(is), fonax_arr(is), &
        bdiag_arr(is)%bunching, bdiag_arr(is)%bunching_phase, bdiag_arr(is)%mean_energy, &
        bdiag_arr(is)%sigma_energy, bdiag_arr(is)%sigma_x, bdiag_arr(is)%sigma_y, &
        bdiag_arr(is)%current, bdiag_arr(is)%n_eff
enddo
call fel_toc (fel_t_stats$)

end subroutine write_diag_rows

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine write_ledger_row ()
!
! Routine to write one energy-ledger row (unaveraged mode). The columns:
!   E_beam      Beam energy relative to the reference, sum(w*(gamma-gamma0))*me
!                 [C*eV = J]. Relative, so the per-record change is not differenced off a
!                 large baseline at its summation-rounding floor (fel-physics.md
!                 sec-numerics).
!   U_field     Total window field energy, sum(P_is)*slice_spacing/c [J].
!   dE_kick     Kick-side change dE_step returned by fel_unavg_step.
!   U_escaped   Cumulative energy slippage carried out of the window, banked at the zero
!                 fill in fel_apply_slippage.
!   U_spont     Cumulative spontaneous deposit energy sum|dE_src|^2. This is the one
!                 field-energy term the kick/deposit duality does not charge to the beam.
!   E_radiated  Cumulative energy the beam radiated under bmad_com's radiation switches.
!                 These are the drawn sums and not expectations, so closure stays exact.
!                 Zero with the switches off.
!
! In a time-dependent run the window is an open system. The exactly closing quantity is
! E_beam + U_field + U_escaped - U_spont + E_radiated. Wakes would be a second, unbanked
! exit channel, and this ledger exists only where they are refused.
!-

subroutine write_ledger_row ()

real(rp) e_beam, u_field, p_mc_l, g0_l
integer is, ip

! The field power per slice comes from the stats loop's fel_field_diag evaluation.
! take_stats_record must have run for this record first.

call fel_tic (fel_t_stats$)
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
call fel_toc (fel_t_stats$)

end subroutine write_ledger_row

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine apply_radiation ()
!
! Routine to apply spontaneous radiation inside FEL elements.
! The routine follows Bmad's global switches bmad_com%radiation_damping_on and
! %radiation_fluctuations_on, as every Bmad tracking path does. Interludes get theirs
! through track1. This covers the custom-tracked FEL step, whose Newton-Lorentz
! equations of motion carry no radiation reaction. Per record:
!
!   damping:      dgamma_j = -(2/3) r_e gamma_j^2 ku^2 aw^2 * INT g^2 ds
!                 Each particle uses its own gamma. The envelope integral runs over the
!                 record, so the unaveraged ramps radiate by their actual strength.
!                 Averaged mode has g = 1.
!   fluctuations: The standard undulator quantum-diffusion variance
!                   1.015e-27 * ku^3 aw^2 F(aw) gamma0^4 per meter
!                 with Genesis's Incoherent.cpp F(aw) fits (the Saldin closed form),
!                 drawn gaussian. Genesis draws uniform scaled by sqrt(3) to reach the
!                 same variance. The sqrt(3) normalizes the uniform draw and must not
!                 appear with a gaussian. One draw per beamlet, as in Genesis, since
!                 independent per-particle kicks would break the quiet start's
!                 per-beamlet harmonic cancellation. Draws are serial, in fixed slice
!                 order, from the one seeded stream. So results do not depend on the
!                 thread count. Only the application runs in parallel.
!
! The drawn energy accumulates into e_rad_cum (the ledger's E_radiated column).
! These are drawn sums and not expectations, so the time-dependent closure stays exact.
!
! Note: In the unaveraged mode this double-counts the grid-captured band, about 3% at
! the reference grid (docs/validation.md), and check_spontaneous.py bounds it. That is
! smaller than the uncertainty of the angular-distribution estimate.
!-

subroutine apply_radiation ()

real(rp) c_loss, fform, sig, gam, p_mc_r, dgam, intg2, s0, g_env, gp_env
integer isl, ipr, ibl, nb, nbl_max
!

if (.not. (bmad_com%radiation_damping_on .or. bmad_com%radiation_fluctuations_on)) return
call fel_tic (fel_t_radiation$)

! The envelope integral over this record. The unaveraged mode takes the midpoint sum on
! its own substep grid. The averaged mode takes dz exactly.

if (fel_mode(ie) == unaveraged$) then
  intg2 = 0
  do ipr = 1, run%ustate%nsub
    s0 = (run%ustate%s - und%dz) + (ipr - 0.5_rp) * run%ustate%dsub
    g_env = fel_und_envelope(run%ustate%l, run%ustate%l_ramp, s0, gp_env)
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

nb = max(1, fbeam%beamlet_size)
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
call fel_toc (fel_t_radiation$)

end subroutine apply_radiation

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine progress_line (at_element_end, i_step, n_step)
!
! Routine to print one row of the progress table: where the walk is and what the light
! and beam are doing. The slow modes (the unaveraged mode costs about 30x the averaged,
! docs/user-guide.md) then show progress. Element boundaries always print. Inside an
! element a 2 s wall clock throttle keeps fast runs quiet.
!
! The row is for a human (docs/user-guide.md). Values are SI-prefixed so a column of
! them lines up. The element name is last so every numeric column stays fixed. The
! one-step elements get no step field. The files carry full precision, and nothing
! parses this row.
!
! The numbers come from the stats record just taken. When the comb keeps no records
! they come from the element-end row, which is filled before this is called. Between
! element ends in that mode the last element-end values are carried forward.
!-

subroutine progress_line (at_element_end, i_step, n_step)

logical at_element_end
integer i_step, n_step
integer(8) now
real(rp) elapsed, pow, ene, bun
character(300) line
character(12) step_str, ela_str

!

call system_clock (now)
if (.not. at_element_end .and. real(now - prog_count_last, rp) / prog_rate < 2.0_rp) return
prog_count_last = now
elapsed = real(now - prog_count0, rp) / prog_rate

if (stats%irec > 0) call fel_stats_exit_light (stats, prog_power, prog_energy, prog_bunch)
pow = prog_power;  ene = prog_energy;  bun = prog_bunch

! One-step elements (every interlude) say nothing with "1/1".

step_str = ''
if (n_step > 1) write (step_str, '(i4, a, i0)') i_step, '/', n_step

if (elapsed < 3600) then
  write (ela_str, '(i0, a, i2.2)') int(elapsed) / 60, ':', mod(int(elapsed), 60)
else
  write (ela_str, '(i0, a, i2.2, a, i2.2)') int(elapsed) / 3600, ':', &
        mod(int(elapsed) / 60, 60), ':', mod(int(elapsed), 60)
endif

! Header once, then the rows. s is the arc length along the reference orbit (meters).
! Through a bending break it exceeds the chord the light takes, and light_corr accounts
! for the difference.

if (.not. prog_header_done) then
  call out_io (s_blank$, r_name, &
        '     %        s       ele       step          power         energy    <|b|>  elapsed  element')
  prog_header_done = .true.
endif

write (line, '(f7.1, f10.3, i7, a, i0, 2x, a9, 1x, a14, 1x, a14, f9.4, a8, 2x, a)') &
      100 * z_now / lat_length, z_now, ie, '/', branch%n_ele_track, adjustr(step_str(1:9)), &
      fel_si_str(pow, 'W'), fel_si_str(ene, 'J'), bun, trim(ela_str), trim(ele%name)
call out_io (s_blank$, r_name, trim(line))

end subroutine progress_line

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine take_stats_record (with_angles)
!
! Routine to take one stats record at the current z_now through fel_stats_record.
! with_angles fills the field theta moments (element ends). An error sets the host's
! err_flag.
!-

subroutine take_stats_record (with_angles)

logical with_angles, serr

call fel_tic (fel_t_stats$)
call fel_stats_record (stats, fbeam, ffield, z_now, ie, with_angles, bdiag_arr, fpow_arr, fonax_arr, serr)
call fel_toc (fel_t_stats$)
if (serr) then
  err_flag = .true.;  return
endif

! The frame series follows the comb. A frame carries the index of the row just taken,
! so the two read on one axis. Every caller of this routine is a comb position. On the
! device the state came back to the host before the row was taken.

if (run%global%dump_at_comb) then
  call fel_tic (fel_t_dumps$)
  call fel_dump_frame (run, max(ie, 0), with_angles, serr)
  call fel_toc (fel_t_dumps$)
  if (serr) then
    err_flag = .true.;  return
  endif
endif
end subroutine take_stats_record

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine end_of_element ()
!
! Routine to handle the end of an element: the evaluated bunch_params (the Tao
! end-of-element pattern) and any requested dumps.
! Mid-run field dumps unrotate exactly as the final dump does. The rotation is a gauge,
! and re-zeroing slip%first keeps the mapping consistent.
!-

subroutine end_of_element ()

logical eerr
character(500) fname

!

call fel_tic (fel_t_ele_end$)
call fel_stats_element_end (stats, fbeam, z_now, eerr)
call fel_toc (fel_t_ele_end$)
if (eerr) then
  err_flag = .true.;  return
endif
if (dump_beam_here(ie) .or. dump_field_here(ie)) call fel_tic (fel_t_dumps$)
if (dump_beam_here(ie)) then
  write (fname, '(2a, i0, 2a)') trim(out_root), '-at', ie, '-', trim(ele%name)
  call fel_dump_beam (run, ele, trim(fname), eerr, s_pos = ele%s)
  if (eerr) then
    err_flag = .true.;  return
  endif
endif

if (dump_field_here(ie)) then
  write (fname, '(2a, i0, 2a)') trim(out_root), '-at', ie, '-', trim(ele%name)
  call fel_dump_field_set (run, trim(fname), eerr, ix_done = ie)
  if (eerr) then
    err_flag = .true.;  return
  endif
endif
if (dump_beam_here(ie) .or. dump_field_here(ie)) call fel_toc (fel_t_dumps$)

end subroutine end_of_element

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine apply_slippage_banked (slippage)
!
! Routine to apply slippage to every field of the set. The escape bank is drained when
! the escaped field is kept.
! Every field slips by the same amount in fundamental wavelengths (one window, lockstep
! rotation). The harm argument only sets the escape bank's slice light-time for a
! harmonic field.
!-

subroutine apply_slippage_banked (slippage)

real(rp) slippage
integer ihh
logical err_b

! While the field is resident the rotation runs against the device records. The ring
! bookkeeping is fel_apply_slippage's, per member as below. The transmitted slice's
! energy comes back at that member's own wavelength and its device planes are zeroed.
! The escape bank is refused with the device, so nothing is drained.

if (run%dev%resident) then
  call fel_tic (fel_t_slippage$)
  do ihh = 1, n_harm
    call fel_device_apply_slippage (run%dev, ihh, ffield(ihh)%slip, ffield(ihh)%wf, slippage, &
                                    harm = ffield(ihh)%harm)
  enddo
  call fel_toc (fel_t_slippage$)
  return
endif

call fel_tic (fel_t_slippage$)
if (keep_escaped_field) then
  do ihh = 1, n_harm
    call fel_apply_slippage (ffield(ihh)%slip, ffield(ihh)%wf, slippage, ffield(ihh)%bank, &
                             ffield(ihh)%harm, err_b)
    if (err_b) then
      err_flag = .true.;  return
    endif
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
call fel_toc (fel_t_slippage$)

end subroutine apply_slippage_banked

!------------------------------------------------------------------------------
!+
! Subroutine set_field_position ()
!
! Routine to set every field's s_position to z_now.
! The field's beamline position is the lattice's, which the walk owns. wavefront_drift
! advances s_position by the optical path it propagated. Through a geometry break that
! path is the chord, while the lattice position advances by the arc.
!-

subroutine set_field_position ()

integer ihh

!

do ihh = 1, n_harm
  ffield(ihh)%wf%s_position = z_now
enddo

end subroutine set_field_position

end subroutine track_fel_line

end module fel_track_line_mod
