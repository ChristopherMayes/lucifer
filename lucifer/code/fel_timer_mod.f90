!+
! Module fel_timer_mod
!
! Wall clock time of a run, accumulated per phase and printed in the footer. The phases
! partition the walk. Every timed region belongs to exactly one phase, so the fractions
! sum to the walk. The leftover gets its own row, named unaccounted.
!
! The timers sit at region boundaries and never inside a parallel loop. A clock call
! inside the slice loop would time the loop it perturbs. The phases that interleave
! per slice (the deposit against its FFT, the field gather against the RK4) come from a
! sampling profiler instead. doc/performance.md gives that split beside this one.
!
! The timers are always on. Each region costs two system_clock calls against
! milliseconds of work, which is below the run-to-run spread in doc/performance.md.
!
! State is module level, like the FFT plan cache and the field-solve kernel cache. The
! phases accumulate from fel_timer_reset, so a driver walking the lattice twice sees the
! total of both walks. Nothing here is read by anything but the report, so the numbers
! never enter the physics. lucifer_smoke_test checks this: two walks in one process
! must agree bit for bit with two processes.
!
! Note: An embedding program that never calls fel_timer_reset gets no table and no
! error. The clock rate is then zero and fel_timer_write returns at once. An error
! return can leave a phase open. That is harmless since a failed run prints no footer.
!-

module fel_timer_mod

use precision_def
use output_mod

implicit none

! The phases, in print order and grouped by where they act. fel_t_walk$ = 0 is the whole
! walk and is the denominator. It is not a leaf.

integer, parameter :: fel_t_walk$ = 0

integer, parameter :: fel_t_und_prep$ = 1        ! Per-step serial prep: plans, kernel cache.
integer, parameter :: fel_t_sc_profile$ = 2      ! The long-range space-charge profile.
integer, parameter :: fel_t_particles$ = 3       ! Transverse maps, RK4, per-slice wake.
integer, parameter :: fel_t_source$ = 4          ! The coherent source prep and kappa fit.
integer, parameter :: fel_t_solve$ = 5           ! Deposit and field solve, per slice.
integer, parameter :: fel_t_unavg$ = 6           ! The unaveraged step, whole.

integer, parameter :: fel_t_slippage$ = 7        ! The slippage rotation and the escape bank.
integer, parameter :: fel_t_radiation$ = 8       ! Spontaneous radiation inside FEL elements.
integer, parameter :: fel_t_migration$ = 9       ! Slice migration and its wake recompute.
integer, parameter :: fel_t_wake$ = 10           ! Whole-window wake kicks.

integer, parameter :: fel_t_seam$ = 11           ! The Bmad-model interlude's slice tracking.
integer, parameter :: fel_t_drift$ = 12          ! The field's drift through a break.
integer, parameter :: fel_t_gen_int$ = 13        ! The transcribed Genesis interlude.

integer, parameter :: fel_t_stats$ = 14          ! The stats record, the diag and ledger rows.
integer, parameter :: fel_t_ele_end$ = 15        ! The element-end bunch parameters.
integer, parameter :: fel_t_dumps$ = 16          ! Mid-run beam and field dumps.

integer, parameter :: fel_t_device$ = 17         ! The resident device step, whole.

integer, parameter :: fel_t_n$ = 17

character(20), parameter :: fel_phase_name(0:fel_t_n$) = [character(20) :: &
      'walk total', &
      'undulator prep', 'space charge', 'particle push', 'coherent source', &
      'field solve', 'unaveraged step', &
      'slippage', 'radiation', 'migration', 'wake', &
      'seam interlude', 'field drift', 'genesis interlude', &
      'stats and diag', 'element end', 'dumps', 'device step']

integer(8), private :: t_open(0:fel_t_n$) = 0
integer(8), private :: t_acc(0:fel_t_n$) = 0
integer(8), private :: n_region(0:fel_t_n$) = 0
integer(8), private :: t_run0 = 0
integer(8), private :: t_rate = 0

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_timer_reset ()
!
! Routine to zero every phase and start the run clock.
! The driver calls this once before it reads its input. So the run total covers the
! whole process, and the work outside the walk is the difference.
!-

subroutine fel_timer_reset ()

t_open = 0
t_acc = 0
n_region = 0
call system_clock (t_run0, t_rate)

end subroutine fel_timer_reset

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_tic (iphase)
!
! Routine to open a phase.
!
! Note: Call this outside every parallel region. A phase must be closed by fel_toc on
! every path that reaches the footer.
!
! Input:
!   iphase -- integer: One of the fel_t_...$ phase indices.
!-

subroutine fel_tic (iphase)

integer iphase

call system_clock (t_open(iphase))

end subroutine fel_tic

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_toc (iphase)
!
! Routine to close a phase. The elapsed clock is added to the phase total and the region
! is counted.
!
! Input:
!   iphase -- integer: The phase fel_tic opened.
!-

subroutine fel_toc (iphase)

integer iphase
integer(8) now

call system_clock (now)
t_acc(iphase) = t_acc(iphase) + (now - t_open(iphase))
n_region(iphase) = n_region(iphase) + 1

end subroutine fel_toc

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_timer_seconds (iphase) result (secs)
!
! Routine to give one phase's accumulated seconds. fel_t_walk$ gives the walk.
!
! Input:
!   iphase -- integer: One of the fel_t_...$ phase indices.
!
! Output:
!   secs   -- real(rp): Seconds accumulated in that phase.
!-

function fel_timer_seconds (iphase) result (secs)

integer iphase
real(rp) secs

secs = 0
if (t_rate > 0) secs = real(t_acc(iphase), rp) / t_rate

end function fel_timer_seconds

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_timer_write (r_name)
!
! Routine to print the phase table.
! Each phase this run entered gets its region count, its seconds and its share of the
! walk. The unaccounted remainder, the walk and the run follow. Phases the run never
! entered are left out.
!
! The phases partition the walk, so the unaccounted row is the walk's own overhead: the
! element loop, the comb tests, the wake and space-charge setup per element, and the
! timers themselves. A large remainder means a phase is missing a timer.
!
! The "= FEL step" row is derived. It is the sum of the phases inside the step, which
! is the share the undulator segments own. It is left out of the leaf sum, so the leaf
! rows still add to the walk minus the remainder.
!
! Input:
!   r_name -- character(*): The caller's name, for out_io.
!-

subroutine fel_timer_write (r_name)

character(*) r_name
character(200) line
character(20) label
real(rp) walk, run, leaf_sum, frac, step
integer(8) now
integer ip

!

if (t_rate == 0) return

walk = fel_timer_seconds(fel_t_walk$)
call system_clock (now)
run = real(now - t_run0, rp) / t_rate

leaf_sum = 0
do ip = 1, fel_t_n$
  leaf_sum = leaf_sum + fel_timer_seconds(ip)
enddo

call out_io (s_blank$, r_name, ' Timing      phase                  regions      seconds     of walk')

! The step's own phases, then their subtotal, then everything outside the step. The
! subtotal follows the group because which phases a run enters depends on its mode.

do ip = fel_t_und_prep$, fel_t_unavg$
  call write_phase_row (ip)
enddo

step = 0
do ip = fel_t_und_prep$, fel_t_unavg$
  step = step + fel_timer_seconds(ip)
enddo
if (step > 0) then
  frac = 0
  if (walk > 0) frac = 100 * step / walk
  label = '= FEL step'
  write (line, '(a, a20, a13, f13.3, f11.1, a)') '             ', label, '', step, frac, '%'
  call out_io (s_blank$, r_name, trim(line))
endif

do ip = fel_t_slippage$, fel_t_n$
  call write_phase_row (ip)
enddo

frac = 0
if (walk > 0) frac = 100 * (walk - leaf_sum) / walk
label = 'unaccounted'
write (line, '(a, a20, a13, f13.3, f11.1, a)') '             ', label, '', &
      walk - leaf_sum, frac, '%'
call out_io (s_blank$, r_name, trim(line))

label = 'walk'
write (line, '(a, a20, i13, f13.3, f11.1, a)') '             ', label, &
      n_region(fel_t_walk$), walk, 100.0_rp, '%'
call out_io (s_blank$, r_name, trim(line))

! The run minus the walk is parsing, the beam and field build, the slippage schedule,
! and the final dumps and stats file.

label = 'run'
write (line, '(a, a20, a13, f13.3)') '             ', label, '', run
call out_io (s_blank$, r_name, trim(line))

!------------------------------------------------------------------------------
contains

!+
! Subroutine write_phase_row (iphase)
!
! Routine to print one phase row, or nothing when the run never entered the phase.
! For example, an averaged run prints no unaveraged step row.
!-

subroutine write_phase_row (iphase)

integer iphase
real(rp) f

if (n_region(iphase) == 0) return
f = 0
if (walk > 0) f = 100 * fel_timer_seconds(iphase) / walk
write (line, '(a, a20, i13, f13.3, f11.1, a)') '             ', fel_phase_name(iphase), &
      n_region(iphase), fel_timer_seconds(iphase), f, '%'
call out_io (s_blank$, r_name, trim(line))

end subroutine write_phase_row

end subroutine fel_timer_write

end module fel_timer_mod
