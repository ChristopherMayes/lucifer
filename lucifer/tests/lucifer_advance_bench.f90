!+
! Program lucifer_advance_bench
!
! The FEL-step particle-path microbenchmark: time fel_advance on a synthetic single
! slice at fixed state and decompose the per-particle-step cost into
! the full path, the RK4+ODE alone (no gather, fixed rpart), and the bare sin/cos
! pair. Each optimization lever then gets a number, not a guess. Serial on purpose:
! this measures the per-particle constant. Parallel scaling is the perf harness's job.
!
!   lucifer_advance_bench [npart] [nstep]      (defaults 8192, 400)
!
! Prints ns per particle-step for: full fel_advance, fel_runge_kutta only, sincos
! only, and a checksum (so the compiler cannot elide the work). Genesis parity note:
! the RK4/ODE arithmetic is BeamSolver::RungeKutta/ODE verbatim, so any cost ratio to
! a clang-built twin of the same loop is compiler/libm, not physics.
!-

program lucifer_advance_bench

use fel_track_mod

implicit none

type (fel_und_struct) und
type (fel_beam_struct), target :: beam
type (fel_field_struct), target :: ff(1)
type (fel_collective_struct) coll
type (fel_slice_struct), pointer :: sl

integer :: npart = 8192, nstep = 400
integer(8) c0, c1, crate
integer ip, istep, ng, n_arg
real(rp) t_full, t_rk, t_sc, delz, phi0, chk
real(rp) gamma, theta, btpar, ez_ip, s, c
complex(rp) rpart
character(20) arg

!

n_arg = command_argument_count()
if (n_arg >= 1) then
  call get_command_argument (1, arg);  read (arg, *) npart
endif
if (n_arg >= 2) then
  call get_command_argument (2, arg);  read (arg, *) nstep
endif

! A synthetic but physically sane state: the benchmark tiers' helical segment
! (gamma0 11357.82, lambda 1e-10, aw 0.84853, lambda_u 0.015), one slice, a seeded
! Gaussian-ish field so the gather reads real numbers.

und%ku = twopi / 0.015_rp
und%aw = 0.84853_rp
und%helical = .true.
und%kx = 0.5_rp * und%ku**2
und%ky = 0.5_rp * und%ku**2
und%pol = [cmplx(1.0_rp, 0.0_rp, rp), cmplx(0.0_rp, -1.0_rp, rp)] / sqrt(2.0_rp)
und%nstep = 1
und%dz = 0.015_rp
delz = 0.015_rp

beam%p0c = sqrt(11357.82_rp**2 - 1) * m_electron
beam%phi0 = 0
beam%wavelength = 1e-10_rp
beam%slice_spacing = 1e-10_rp
beam%beamlet_size = 8
allocate (beam%slice(1))
sl => beam%slice(1)
call fel_slice_reallocate (sl, npart)
sl%n = npart

call ran_seed_put (777)
do ip = 1, npart
  call ran_gauss (sl%x(ip));   sl%x(ip) = 3e-5_rp * sl%x(ip)
  call ran_gauss (sl%y(ip));   sl%y(ip) = 3e-5_rp * sl%y(ip)
  call ran_gauss (sl%px(ip));  sl%px(ip) = 1e-6_rp * sl%px(ip)
  call ran_gauss (sl%py(ip));  sl%py(ip) = 1e-6_rp * sl%py(ip)
  call ran_gauss (sl%pz(ip));  sl%pz(ip) = 1e-4_rp * sl%pz(ip)
  sl%z(ip) = 1e-10_rp * (real(ip, rp) / npart - 0.5_rp)
  sl%weight(ip) = 1e-15_rp
enddo

ng = 151
ff(1)%harm = 1
call wavefront_init (ff(1)%wf, ng, ng, 1, 2 * 2e-4_rp / (ng - 1), 2 * 2e-4_rp / (ng - 1), &
                     1e-10_rp, 1e-10_rp, 'x', 0.0_rp)
ff(1)%wf%Ex = cmplx(3e5_rp, 1e5_rp, rp)     ! Flat nonzero field: the gather does real work.
ff(1)%slip%timerun = .false.
ff(1)%slip%sample = 1

call system_clock (c0, crate)

! ---- 1. The full particle path, as the walk runs it (nf == 1 branch).

phi0 = 0.1_rp
call system_clock (c0)
do istep = 1, nstep
  call fel_advance (und, beam, sl, ff, delz, phi0, coll, 1)
enddo
call system_clock (c1)
t_full = real(c1 - c0, rp) / crate

! ---- 2. RK4 + ODE only: fixed gather result, same per-stage arithmetic.

chk = 0
rpart = cmplx(1e-3_rp, 5e-4_rp, rp)
ez_ip = 0
btpar = 1 + und%aw**2
call system_clock (c0)
do istep = 1, nstep
  do ip = 1, npart
    gamma = 11357.82_rp + 1e-3_rp * ip
    theta = 1e-4_rp * ip
    call fel_runge_kutta (delz, twopi / 1e-10_rp, und%ku, btpar, rpart, ez_ip, gamma, theta)
    chk = chk + gamma
  enddo
enddo
call system_clock (c1)
t_rk = real(c1 - c0, rp) / crate

! ---- 3. The bare transcendental pair, 4 per particle-step as the ODE calls them.

call system_clock (c0)
do istep = 1, nstep
  do ip = 1, npart
    theta = 1e-4_rp * ip + 1e-7_rp * istep
    s = sin(theta);  c = cos(theta)
    chk = chk + s + c
    s = sin(theta + 0.1_rp);  c = cos(theta + 0.1_rp)
    chk = chk + s + c
    s = sin(theta + 0.2_rp);  c = cos(theta + 0.2_rp)
    chk = chk + s + c
    s = sin(theta + 0.3_rp);  c = cos(theta + 0.3_rp)
    chk = chk + s + c
  enddo
enddo
call system_clock (c1)
t_sc = real(c1 - c0, rp) / crate

print '(a, i0, a, i0, a)', 'lucifer_advance_bench: ', npart, ' particles x ', nstep, ' steps, serial'
print '(a, f8.2, a)', '  full fel_advance:      ', 1e9_rp * t_full / npart / nstep, ' ns/particle-step'
print '(a, f8.2, a)', '  fel_runge_kutta only:  ', 1e9_rp * t_rk / npart / nstep, ' ns/particle-step'
print '(a, f8.2, a)', '  4x (sin+cos) only:     ', 1e9_rp * t_sc / npart / nstep, ' ns/particle-step'
print '(a, es12.4)', '  checksum (do not optimize away): ', chk + sum(sl%pz(1:10))

end program lucifer_advance_bench
