!+
! Module fel_cost_mod
!
! What a run will cost, counted before it tracks and estimated from rates measured on
! one machine (doc/performance.md).
!
! The count is exact. Once the setup pass has run, the slices, the macroparticles, the
! integration steps and the grid points are all known, and their product is the work.
! The estimate is not exact, and the line that prints it names the machine it is for.
!
! The model is two terms, which is what the measurement supports. The particle push is
! linear in particle-steps. The field solve is a transform that does not care how many
! particles deposited into the grid, going as ngrid^2 log2(ngrid) per slice per step,
! with a deposit linear in particle-steps beside it. The 512 and 8192 particle rows of
! the performance page separate the two, since the push moves by a factor of 14 across
! them while the transform moves by 12 percent.
!
! Rates below are per thread. The page's tables were taken at twelve threads, and the
! thread scan there inverts to a serial fraction of 2.8 percent, which reproduces the
! measured speedup at 2, 4, 8 and 12 threads to better than 2 percent. The parallelism
! is over slices, so a window of one slice gets nothing from threads and the effective
! thread count is the smaller of the two.
!
! Against the eight configurations the performance page measures, the estimate lands
! within 2 percent. That is the calibration set and not a claim about another machine
! or another line. On the examples and the comparison tiers the estimate runs low by 20
! to 40 percent, and part of that is dated: the page's runs were taken before the source
! filter became the default, and the filter costs 8 to 14 percent of wall clock
! (doc/fel-physics.md sec-convergence). The estimate is a number to plan a run against,
! so it is checked against the clock in the footer of every run rather than trusted.
!-

module fel_cost_mod

use fel_struct
!$ use omp_lib

implicit none

! The recording machine and the day the rates were fitted, printed with the estimate so
! that a number read on another machine is read as what it is.

character(*), parameter :: fel_cost_machine$ = 'Apple M3 Max, 12 cores, 2026-09-07'

! Fitted from doc/performance.md, the phase profile and the particles-against-the-grid
! tables, production build, divided out to one thread.

real(rp), parameter :: fel_cost_push$ = 1.96e-7_rp     ! s per particle-step.
real(rp), parameter :: fel_cost_deposit$ = 2.09e-8_rp  ! s per particle-step.
real(rp), parameter :: fel_cost_solve$ = 2.54e-9_rp    ! s per grid point per log2(ngrid) per step.

! The thread scan's serial fraction, and the share of the walk the push and the solve
! own. The share is 83.6, 83.5 and 82.9 percent over a factor of 16 in the load and
! 82.9 at 504 slices, so one constant carries it.

real(rp), parameter :: fel_cost_serial$ = 0.028_rp
real(rp), parameter :: fel_cost_walk_share$ = 0.83_rp

! The device against the CPU, from the page's own table: the largest window measured is
! 12.4 times the twelve-thread CPU, and the smallest deck measures the host dispatch
! floor at 0.037 s over 89 steps. The four device points cannot carry a two-term fit of
! their own, since the two largest share a particle count per slice and the smallest is
! the floor, so the device takes the CPU model divided by the measured ratio with the
! floor added. The ratio runs from 7.3 to 12.4 over the cases measured, so this term is
! the weakest one here and the estimate says so by naming the backend.

real(rp), parameter :: fel_cost_device$ = 12.4_rp
real(rp), parameter :: fel_cost_dispatch$ = 4.2e-4_rp  ! s per step, host side.
integer, parameter :: fel_cost_device_threads$ = 12    ! The CPU thread count that ratio is against.

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_cost_count (run, nstep, n_particle, n_pstep, n_gpt, ngrid)
!
! Routine to count the work of a run from the state the setup pass built. Exact: every
! number here is known before the first step.
!
! Input:
!   run        -- fel_run_struct: The run, after fel_setup_schedule.
!
! Output:
!   nstep      -- integer: FEL integration steps over the walk, summed over the segments
!                   the walk covers. An unaveraged segment counts its own substeps,
!                   which is what it integrates and deposits on.
!   n_particle -- integer: Macroparticles in the window.
!   n_pstep    -- real(rp): Particle-steps, the product of the two. Real because it
!                   passes 2^31 on an ordinary time-dependent run.
!   n_gpt      -- real(rp): Transverse grid points the field solve touches per step,
!                   over the window and over every member and plane of the field set.
!   ngrid      -- integer: Transverse points on a side of one grid.
!   nstep_unavg -- integer: How many of nstep are unaveraged substeps.
!-

subroutine fel_cost_count (run, nstep, n_particle, n_pstep, n_gpt, ngrid, nstep_unavg)

type (fel_run_struct), target :: run
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
real(rp) n_pstep, n_gpt, lambda_w
integer nstep, n_particle, ngrid, nstep_unavg, ie, n_ele, nsub, n_grid(3), npol

!

! Only the elements the walk covers. global%track_start and track_end bound it, and a
! deck that stops after the second segment pays for two segments and not for the line.

branch => run%lat%branch(0)
nstep = 0;  nstep_unavg = 0
do ie = run%i_start, run%i_end
  if (.not. run%is_fel(ie)) cycle
  ele => branch%ele(ie)
  n_ele = max(1, nint(ele%value(num_steps$)))

  ! An unaveraged segment resolves the quiver, so its integration step is a fraction of
  ! a period rather than the record step, and it pushes and deposits on every substep.
  ! fel_unavg_setup builds the same count from the same two numbers.

  if (run%fel_mode(ie) == unaveraged$) then
    lambda_w = ele%value(l_period$)
    nsub = 1
    if (lambda_w > 0) nsub = max(1, nint((ele%value(l$) / n_ele) * &
                                         run%global%unaveraged_steps_per_period / lambda_w))
    nstep_unavg = nstep_unavg + n_ele * nsub
    nstep = nstep + n_ele * nsub
  else
    nstep = nstep + n_ele
  endif
enddo

n_particle = sum(run%fbeam%slice%n)
n_pstep = real(n_particle, rp) * nstep

n_grid = wavefront_shape(run%ffield(1)%wf)
ngrid = n_grid(1)
npol = 1
if (run%two_pol) npol = 2
n_gpt = real(run%nslice, rp) * n_grid(1) * n_grid(2) * run%n_harm * npol

end subroutine fel_cost_count

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_cost_estimate (run, secs, backend)
!
! Routine to estimate the wall clock of a run's walk on the recording machine.
!
! The estimate covers the walk, which is the tracking. The parse, the beam and field
! build and the final dumps are outside it and are 1 percent of the walk on the page's
! own runs.
!
! Input:
!   run     -- fel_run_struct: The run, after fel_setup_schedule.
!
! Output:
!   secs    -- real(rp): Estimated walk seconds.
!   backend -- character(*): 'CPU' or the device name, whichever the estimate is for.
!-

subroutine fel_cost_estimate (run, secs, backend)

type (fel_run_struct), target :: run
real(rp) secs, n_pstep, n_gpt, t_one, speedup, f_unavg
integer nstep, n_particle, ngrid, nstep_unavg, n_omp, n_eff
character(*) backend

!

call fel_cost_count (run, nstep, n_particle, n_pstep, n_gpt, ngrid, nstep_unavg)

! One thread, the two terms. Every member and plane of the field set is already counted
! into n_gpt, and each carries the same log of its own grid.

t_one = (fel_cost_push$ + fel_cost_deposit$) * n_pstep
if (ngrid > 1) t_one = t_one + fel_cost_solve$ * n_gpt * log(real(ngrid, rp)) / log(2.0_rp) * nstep

n_omp = 1
!$ n_omp = omp_get_max_threads()
n_eff = max(1, min(n_omp, run%nslice))
speedup = 1 / (fel_cost_serial$ + (1 - fel_cost_serial$) / n_eff)

if (run%dev%on) then
  backend = trim(run%global%device)
  speedup = 1 / (fel_cost_serial$ + (1 - fel_cost_serial$) / fel_cost_device_threads$)
  secs = t_one / speedup / fel_cost_device$ / fel_cost_walk_share$ + fel_cost_dispatch$ * nstep
else
  backend = 'CPU'
  secs = t_one / speedup / fel_cost_walk_share$
endif

end subroutine fel_cost_estimate

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_cost_secs_str (secs) result (str)
!
! Routine to write a number of seconds for a human. Seconds are the unit either way, so
! fel_si_str is not used here: its bare-unit row leaves two spaces before the unit, and
! its prefixes would report a short run in milliseconds where the line beside it is in
! seconds.
!
! Input:
!   secs -- real(rp): Seconds.
!
! Output:
!   str  -- character(16): The number and its unit, left justified.
!-

function fel_cost_secs_str (secs) result (str)

real(rp) secs
character(16) str

!

if (secs >= 1) then
  write (str, '(f0.1, a)') secs, ' s'
else
  write (str, '(f6.3, a)') secs, ' s'    ! f0.3 drops the leading zero.
  str = adjustl(str)
endif

end function fel_cost_secs_str

end module fel_cost_mod
