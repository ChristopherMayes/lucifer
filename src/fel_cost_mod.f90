!+
! Module fel_cost_mod
!
! Work count and wall clock estimate of a run, made before it tracks. The rates come
! from one machine (docs/performance.md).
!
! The count is exact. After the setup pass the slices, the macroparticles, the
! integration steps and the grid points are all known. Their product is the work.
! The estimate is not exact, and the line that prints it names the machine it is for.
!
! The model has two terms. The particle push is linear in particle-steps. The field
! solve is a transform that goes as ngrid^2 log2(ngrid) per slice per step, whatever the
! particle count. A deposit linear in particle-steps sits beside it. The 512 and 8192
! particle rows of docs/performance.md separate the two terms: the push moves by a factor
! of 14 between them and the transform by 8 percent.
!
! The rates below are per thread. The page's tables were taken at twelve threads. The
! thread scan there gives a serial fraction of 2.1 percent. The parallelism is over
! slices, so the effective thread count is the smaller of the thread and slice counts.
! A window of one slice gets nothing from threads.
!
! tests/scripts/fit_cost_rates.py fits the rates from the logs of the runs the page
! records: sixteen runs, three particle counts and two slice counts at five thread
! counts. The worst residual on the walk is 16.5 percent over that set, on the thread
! sweep's four- and eight-thread rows, where one serial fraction cannot carry a sweep
! whose efficiency is 92 percent at four threads and 94 at twelve. The footer of every
! run prints the walk beside the estimate.
!
! Note: The rates hold for the recording machine and line. Another machine or another
! line can differ.
!-

module fel_cost_mod

use fel_struct
!$ use omp_lib

implicit none

! The recording machine and the day the rates were fitted, printed with the estimate.

character(*), parameter :: fel_cost_machine$ = 'Apple M3 Max, 12 cores, 2026-10-24'

! Fitted from the phase profile and the particles-against-the-grid tables of
! docs/performance.md (production build, gfortran 15.3), divided out to one thread. The
! source filter is on in every run of those tables since it is the default. The filter
! adds a forward transform of the source per field per step, which the measurement of
! 2026-09-07 put at 58 percent of the transform rate, with the push and the deposit
! left alone.

real(rp), parameter :: fel_cost_push$ = 2.03e-7_rp     ! s per particle-step.
real(rp), parameter :: fel_cost_deposit$ = 2.08e-8_rp  ! s per particle-step.
real(rp), parameter :: fel_cost_solve$ = 4.68e-9_rp    ! s per grid point per log2(ngrid) per step.

! The thread scan's serial fraction, and the share of the walk the push and the solve
! own. The share is 88.5, 87.5 and 85.9 percent over a factor of 16 in particles per
! slice and 88.1 percent at 504 slices, so one constant carries it.

real(rp), parameter :: fel_cost_serial$ = 0.0105_rp
real(rp), parameter :: fel_cost_walk_share$ = 0.883_rp

! The device against the CPU, from the device table of docs/performance.md. The largest
! window there is 9.1 times the twelve-thread CPU. The smallest deck gives the host
! dispatch floor: 0.045 s over 89 steps.
! The four device points cannot carry a two-term fit of their own. The two largest share
! a particle count per slice and the smallest is the floor. So the device takes the CPU
! model divided by the ratio, plus the floor.
!
! Note: The ratio runs from 6.1 to 9.1 over the device table. This is the least certain
! term, and the printed estimate names the backend.

real(rp), parameter :: fel_cost_device$ = 9.1_rp
real(rp), parameter :: fel_cost_dispatch$ = 5.1e-4_rp  ! s per step, host side.
integer, parameter :: fel_cost_device_threads$ = 12    ! The CPU thread count that ratio is against.

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_cost_count (run, nstep, n_particle, n_pstep, n_gpt, ngrid, nstep_unavg)
!
! Routine to count the work of a run from the state the setup pass built.
! The count is exact since every number is known before the first step.
!
! Input:
!   run         -- fel_run_struct: The run, after fel_setup_schedule.
!
! Output:
!   nstep       -- integer: FEL integration steps, summed over the segments the walk
!                    covers. An unaveraged segment counts its substeps, since it
!                    integrates and deposits on each.
!   n_particle  -- integer: Macroparticles in the window.
!   n_pstep     -- real(rp): Particle-steps, the product of the two. Real since it
!                    passes 2^31 on an ordinary time-dependent run.
!   n_gpt       -- real(rp): Transverse grid points the field solve touches per step,
!                    over the window and over every member and plane of the field set.
!   ngrid       -- integer: Transverse points on a side of one grid.
!   nstep_unavg -- integer: How many of nstep are unaveraged substeps.
!-

subroutine fel_cost_count (run, nstep, n_particle, n_pstep, n_gpt, ngrid, nstep_unavg)

type (fel_run_struct), target :: run
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
real(rp) n_pstep, n_gpt, lambda_w
integer nstep, n_particle, ngrid, nstep_unavg, ie, n_ele, nsub, n_grid(3), npol

!

! Only the elements the walk covers, as bounded by global%track_start and track_end.
! For example, a deck that stops after the second segment counts two segments.

branch => run%lat%branch(0)
nstep = 0;  nstep_unavg = 0
do ie = run%i_start, run%i_end
  if (.not. run%is_fel(ie)) cycle
  ele => branch%ele(ie)
  n_ele = max(1, nint(ele%value(num_steps$)))

  ! An unaveraged segment resolves the quiver. Its integration step is a fraction of a
  ! period, and it pushes and deposits on every substep. fel_unavg_setup builds the same
  ! count from the same two numbers.

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
! build and the final dumps are outside it. They are 1 percent of the walk on the runs
! of docs/performance.md.
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

! One thread, the two terms. Every member and plane of the field set is already in n_gpt,
! and each has the same log2(ngrid).

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
! Routine to write a number of seconds for a human.
! fel_si_str is not used. Its bare-unit row leaves two spaces before the unit. And its
! prefixes would put a short run in milliseconds beside a line in seconds.
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
