!+
! Module fel_fp32_mod
!
! The single-precision particle path and the lockstep instrument that prices it
! (doc/performance.md, doc/validation.md). Apple GPUs have no FP64, so before any
! kernel is written a single-precision form of the averaged FEL advance must exist
! with its divergence from the FP64 path a measured number. Nothing outside the
! instrument reads the FP32 state, and the FP64 path is untouched: this module is a
! measuring device, not a production mode.
!
! The FP32 state is the packed six-vector (x, px, y, py, dzr, pz), 24 bytes per
! particle. dzr is the longitudinal residual of the split z = z_ref + dzr, with z_ref
! one FP64 number per slice. Storing Bmad's z absolutely in FP32 fails silently (the
! per-step increment rounds to zero and bunching never forms), so the residual form is
! load-bearing and the runtime guard below watches it. The reference must also move:
! the beam's common z drift is secular, a static reference lets the residual grow until
! its ulp swallows the physics, and the guard caught exactly that on this instrument's
! first run (step 514 of the steady-state example, 28 ulps against the floor of 32).
! z_ref therefore refreshes to the slice's FP64 mean each step, one host number per
! slice, and in freerun mode the persistent residuals rebase across the reference move
! by fel_fp32_renorm, which is the migration operation doing double duty.
!
! Reformulations, each forced by an FP32 quantum and each below FP64 notice:
!
!   energy   The working variable is goff = gamma - gamma0. The conversion is
!            goff = pz * p0_mc^2/gamma0, exact to O(1/gamma^2) ~ 4e-9, below FP32
!            resolution. Absolute FP32 gamma has a quantum of 1.35e-3 at these
!            energies against a per-step change of 3.0e-4, and sqrt(gamma^2 - 1)
!            loses the 1 entirely (the quantum of gamma^2 is about 15).
!   phase    The working variable is delta = theta - phi0, formed as ks * dzr
!            (the 1/beta corrections are ~4e-9 of theta, below FP32 resolution and
!            part of what the instrument measures). The phase factor is
!            e^{-i theta} = e^{-i phi0} * e^{-i delta}: the base rotator is one FP64
!            evaluation per slice, the per-particle sincos is FP32 on the small
!            angle. This is the angle-addition form a device kernel would use.
!   detuning The FP64 ODE forms sqrt(1 - btper0/gamma^2), which is 1 minus ~1.3e-8:
!            in FP32 that is exactly 1 and the phase equation loses the energy
!            dependence completely. The FP32 ODE forms the detuning as the
!            difference 0.5*ks*(qres - q) with q = btper0/gamma^2 and
!            qres = 2*ku/ks, so the resonant cancellation happens between two
!            small like quantities. The remaining floor is btper0's own FP32
!            representation, and the recorded levels carry it.
!
! The instrument runs both precisions from a shared state at fel_advance's own
! sequence point: the FP64 slice is copied before its
! advance, the twin advances the FP32 image of that copy, and the comparison reads
! both results before anything else touches the slice. In lockstep mode the FP32
! state is rebuilt from FP64 every step, so a wrong formula shows as a jump rather
! than as growth. In freerun mode dzr and pz persist across steps and the compounding
! rate is measured; the transverse coordinates are refreshed from FP64 in both modes,
! since the transverse maps stay FP64 and the twin covers what the FEL exchange
! evolves. The source phasor sum(w * awloc * e^{-i theta}) is compared at the
! same point; the production deposit runs one FP64 transverse half-step later, which
! is precision-neutral between the sides.
!
! Everything here is serial-safe inside the caller's parallel slice loop: the twin
! touches only its slice's state and its slice's row of the divergence arrays, so the
! instrument's own output is identical at any thread count.
!-

module fel_fp32_mod

use fel_beam_mod
use wavefront_mod

implicit none

integer, parameter :: fel_fp32_nq$ = 7   ! x, px, y, py, pz, theta, phasor.

character(8), parameter :: fel_fp32_qname(fel_fp32_nq$) = &
      [character(8) :: 'x', 'px', 'y', 'py', 'pz', 'theta', 'phasor']

!+
! Struct fel_fp32_slice_struct
!
! One slice's packed FP32 particle set: the six-vector plus the working longitudinal
! chart the twin evolves. Allocated to the slice's fill count at first use.
!-

type fel_fp32_slice_struct
  real(sp), allocatable :: x(:), px(:), y(:), py(:)   ! Transverse, Bmad conventions.
  real(sp), allocatable :: dzr(:)                     ! z residual off z_ref [m].
  real(sp), allocatable :: pz(:)                      ! (p - p0)/p0, already an offset.
  real(rp) :: z_ref = 0                               ! The slice's FP64 z reference [m].
end type

!+
! Struct fel_fp32_struct
!
! The instrument's run state: mode, the per-slice FP32 sets, the per-step divergence
! table the parallel loop fills and the serial epilogue reduces, and the guard and
! worst-case accumulators the footer and the check read.
!-

type fel_fp32_struct
  logical :: on = .false.
  logical :: freerun = .false.       ! dzr/pz persist across steps (compounding measured).
  logical :: mutate = .false.        ! Check hook: truncate dzr harder so the check can fail.
  integer :: iu = 0                  ! The .fp32.txt stream.
  integer :: istep = 0               ! Steps instrumented.
  type (fel_fp32_slice_struct), allocatable :: sl32(:)
  real(rp), allocatable :: div_slice(:,:)    ! (nq, nslice): this step, filled in parallel.
  real(rp), allocatable :: ulp_slice(:)      ! (nslice): median per-step |d delta| in ulps.
  real(rp) :: worst(fel_fp32_nq$) = 0        ! Run-level worst per quantity.
  real(rp) :: ulp_min = huge(1.0_rp)         ! Run-level worst (smallest) guard statistic.
  logical :: checked = .false.               ! A configuration check ran (first step).
end type

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_fp32_setup (fp32, mode, mutate, nslice, out_root, err_flag)
!
! Routine to arm the instrument from the input knob: '' or 'off' leaves it dark,
! 'lockstep' rebuilds the FP32 state from FP64 every step, 'freerun' lets the
! longitudinal state carry across steps. Opens the .fp32.txt stream.
!
! Input:
!   mode     -- character(*): global%fp32_check.
!   mutate   -- logical: global%fp32_mutate, the check's own failure hook.
!   nslice   -- integer: Slice count.
!   out_root -- character(*): Output root for the stream file.
!
! Output:
!   fp32     -- fel_fp32_struct: Armed (or left dark).
!   err_flag -- logical: Set True on an unrecognized mode. False otherwise.
!-

subroutine fel_fp32_setup (fp32, mode, mutate, nslice, bucket_shift, out_root, err_flag)

type (fel_fp32_struct) fp32
character(*) mode, out_root
logical mutate, err_flag
integer nslice
real(rp) bucket_shift
real(sp), allocatable :: rres(:), rorig(:)
real(rp) renorm_worst
integer i, nprobe
character(*), parameter :: r_name = 'fel_fp32_setup'

!

err_flag = .false.

select case (mode)
case ('', 'off')
  return
case ('lockstep')
  fp32%freerun = .false.
case ('freerun')
  fp32%freerun = .true.
case default
  call out_io (s_error$, r_name, 'UNRECOGNIZED FP32_CHECK MODE: "' // trim(mode) // '".', &
                                 'THE MODES ARE off, lockstep AND freerun.')
  err_flag = .true.
  return
end select

fp32%on = .true.
fp32%mutate = mutate
allocate (fp32%sl32(nslice))
allocate (fp32%div_slice(fel_fp32_nq$, nslice))
allocate (fp32%ulp_slice(nslice))

open (newunit = fp32%iu, file = trim(out_root) // '.fp32.txt', action = 'write')
write (fp32%iu, '(a)') '# FP32 lockstep instrument: per-step worst relative divergence per quantity,'
write (fp32%iu, '(a)') '# FP32 twin against the FP64 path from a shared state at fel_advance. theta is'
write (fp32%iu, '(a)') '# absolute [rad], guard_ulp is the median per-step residual phase increment in'
write (fp32%iu, '(a)') '# ulps of the residual (the silent-z guard: small means FP32 cannot resolve the'
write (fp32%iu, '(a)') '# step and the run refuses).'
write (fp32%iu, '(a, l1)') '# freerun = ', fp32%freerun

! Migration in the residual representation is a bucket renormalization: a mover's
! residual re-references to the new slice origin by the FP64 bucket shift. The
! round trip (shift out and back) is measured here on residuals spanning the slice,
! in ulps of the residual, and recorded beside the run's own levels.

nprobe = 4097
allocate (rres(nprobe), rorig(nprobe))
do i = 1, nprobe
  rres(i) = real((real(i-1, rp) / (nprobe-1) - 0.5_rp) * bucket_shift, sp)
enddo
rorig = rres
call fel_fp32_renorm (rres, bucket_shift)
call fel_fp32_renorm (rres, -bucket_shift)
renorm_worst = 0
do i = 1, nprobe
  renorm_worst = max(renorm_worst, real(abs(rres(i) - rorig(i)), rp) / &
                     max(real(spacing(max(abs(rorig(i)), real(bucket_shift, sp))), rp), tiny(1.0_rp)))
enddo
write (fp32%iu, '(a, es13.4, a, es13.4, a)') '# renorm_roundtrip_ulp ', renorm_worst, &
      '   (bucket shift ', bucket_shift, ' m)'
write (fp32%iu, '(a)') '#   step          x            px           y            py           pz         theta        phasor      guard_ulp'

end subroutine fel_fp32_setup

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_fp32_renorm (dzr, shift)
!
! Routine to re-reference FP32 residuals across a bucket boundary: the migration
! operation in the split representation. The shift is FP64 (it is the slice bucket,
! exactly beta * slice_spacing) and rounds once into the residual's own precision.
!
! Input:
!   dzr(:) -- real(sp): Residuals off the old reference.
!   shift  -- real(rp): The bucket shift [m], signed.
!
! Output:
!   dzr(:) -- real(sp): Residuals off the new reference.
!-

subroutine fel_fp32_renorm (dzr, shift)

real(sp) dzr(:)
real(rp) shift

dzr = dzr - real(shift, sp)

end subroutine fel_fp32_renorm

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_fp32_twin_slice (fp32, is, x0, y0, px0, py0, z0, pz0, sl, beam, &
!                                    exfld, dx, dy, gridmax, aw, ku, kx, ky, ax, ay, &
!                                    cos_t, sin_t, rtmp, delz, xks, xku, phi0, phi0_new)
!
! Routine to advance one slice's FP32 twin from the shared pre-advance state and fill
! this slice's divergence row. The x0.. arrays are the FP64 state exactly as
! fel_advance received it, sl is the slice after the FP64 advance, and the twin
! mirrors fel_advance's arithmetic under the module header's reformulations:
! fundamental field only, no collective terms (the setup refusals hold that).
!
! Runs inside the caller's parallel loop: everything written is indexed by is.
!-

subroutine fel_fp32_twin_slice (fp32, is, x0, y0, px0, py0, z0, pz0, sl, beam, &
                                exfld, dx, dy, gridmax, aw, ku, kx, ky, ax, ay, &
                                cos_t, sin_t, rtmp, delz, xks, xku, phi0, phi0_new)

type (fel_fp32_struct), target :: fp32
type (fel_fp32_slice_struct), pointer :: t
type (fel_slice_struct) sl
type (fel_beam_struct) beam
integer is
real(rp) x0(:), y0(:), px0(:), py0(:), z0(:), pz0(:)
complex(rp) exfld(:,:)
real(rp) dx, dy, gridmax, aw, ku, kx, ky, ax, ay, cos_t, sin_t, rtmp, delz, xks, xku
real(rp) phi0, phi0_new

! FP32 constants of the step, converted once.

real(sp) aw32, ku32, kx32, ky32, ax32, ay32, ct32, st32, rt32, dz32, ks32, gam032
real(sp) e032, p032, qres32, gmax32, dxg32, dyg32, cret32
complex(sp) base32

! FP32 per-particle working set.

real(sp) x, y, px, py, pz, dzr, goff, delta, gam, awloc, btpar, px_g, py_g
real(sp) wx, wy, q, s_d, c_d, d_new
complex(sp) cpart, rpart, p32sum
complex(rp) p64sum, e_ip
real(rp) wsum
real(rp) gamma0, p0_mc, theta64, beta64, p_mc64, aw64loc, ulp_med, phi_ref, zref_new
real(rp) dstat(fel_fp32_nq$), sc(fel_fp32_nq$), dinc
real(rp), allocatable :: incr(:)
integer ip, n, ix, iy
logical on_grid

!

t => fp32%sl32(is)
n = sl%n
p0_mc = fel_p0_mc(beam)
gamma0 = fel_gamma0(beam)

if (.not. allocated(t%x)) then
  allocate (t%x(n), t%px(n), t%y(n), t%py(n), t%dzr(n), t%pz(n))
  t%z_ref = sum(z0(1:n)) / n
  t%pz = real(pz0(1:n), sp)
  t%dzr = real(z0(1:n) - t%z_ref, sp)
endif

! The shared state: transverse always from FP64 (the transverse maps stay FP64), the
! longitudinal from FP64 in lockstep mode and carried in freerun mode. The FP64
! reference moves to the slice mean each step in both modes: one host number per
! slice, which is the split's working principle. In freerun the persistent residuals
! rebase across the move, the renormalization the migration path uses.

t%x = real(x0(1:n), sp);   t%px = real(px0(1:n), sp)
t%y = real(y0(1:n), sp);   t%py = real(py0(1:n), sp)
if (fp32%freerun) then
  zref_new = sum(z0(1:n)) / n
  call fel_fp32_renorm (t%dzr, zref_new - t%z_ref)
  t%z_ref = zref_new
else
  t%z_ref = sum(z0(1:n)) / n
  t%pz = real(pz0(1:n), sp)
  t%dzr = real(z0(1:n) - t%z_ref, sp)
endif

! The check's own mutation: a deliberately coarser residual, so the recorded level
! must move and the check is falsifiable.

if (fp32%mutate) t%dzr = anint(t%dzr / (256 * spacing(t%dzr))) * (256 * spacing(t%dzr))

! Step constants. e0 = p0^2/gamma0 converts pz to goff = gamma - gamma0, exact to
! O(1/gamma^2) ~ 4e-9 (module header). qres = 2*ku/ks is the resonance in the
! detuning difference. cret = phi0_new - phi0 is the chart-return base, small.

aw32 = real(aw, sp);  ku32 = real(ku, sp);  kx32 = real(kx, sp);  ky32 = real(ky, sp)
ax32 = real(ax, sp);  ay32 = real(ay, sp);  ct32 = real(cos_t, sp);  st32 = real(sin_t, sp)
rt32 = real(rtmp, sp);  dz32 = real(delz, sp);  ks32 = real(xks, sp)
gam032 = real(gamma0, sp);  p032 = real(p0_mc, sp)
e032 = real(p0_mc**2 / gamma0, sp)
qres32 = real(2 * xku / xks, sp)
gmax32 = real(gridmax, sp);  dxg32 = real(dx, sp);  dyg32 = real(dy, sp)
! The chart return is z_ref-free: theta and dzr are both referenced to phi_ref, so the
! reference cancels and cret stays the small per-step phase advance.

cret32 = real(phi0_new - phi0, sp)
phi_ref = phi0 + xks * t%z_ref
base32 = cmplx(cmplx(cos(phi_ref), -sin(phi_ref), rp), kind=sp)   ! e^{-i(phi0 + ks z_ref)}, FP64 once.

allocate (incr(n))
p32sum = 0
p64sum = 0
wsum = 0
dstat = 0

do ip = 1, n

  x = t%x(ip);  y = t%y(ip);  px = t%px(ip);  py = t%py(ip)
  pz = t%pz(ip);  dzr = t%dzr(ip)

  goff = pz * e032
  delta = ks32 * dzr
  gam = gam032 + goff

  awloc = faw32(x, y)
  px_g = px * p032
  py_g = py * p032
  btpar = 1 + px_g*px_g + py_g*py_g + aw32*aw32*awloc*awloc

  call gather32 (x, y, cpart, on_grid)
  if (on_grid) then
    rpart = rt32 * awloc * conjg(cpart)
  else
    rpart = 0
  endif

  call rk32 (goff, delta)

  ! Chart return, the reformulated forms: pz from goff through e0, the residual from
  ! delta with the beta corrections below FP32 resolution (module header).

  d_new = delta
  pz = goff / e032
  dzr = -(cret32 - d_new) / ks32

  incr(ip) = abs(real(d_new - ks32 * t%dzr(ip), rp))
  t%pz(ip) = pz
  t%dzr(ip) = dzr

  ! The FP32 side of the source phasor, at this same sequence point.

  call sincos32 (d_new, s_d, c_d)
  p32sum = p32sum + real(sl%weight(ip), sp) * awloc * base32 * cmplx(c_d, -s_d, sp)

  ! The FP64 comparison values, from the slice fel_advance just advanced.

  p_mc64 = p0_mc * (1 + sl%pz(ip))
  beta64 = p_mc64 / sqrt(p_mc64**2 + 1)
  theta64 = (phi0_new - phi0) + xks * sl%z(ip) / beta64 - xks * t%z_ref   ! delta convention.
  aw64loc = faw64 (sl%x(ip), sl%y(ip))
  e_ip = cmplx(cos(phi_ref + theta64), -sin(phi_ref + theta64), rp)
  p64sum = p64sum + sl%weight(ip) * aw64loc * e_ip
  wsum = wsum + sl%weight(ip) * aw64loc

  ! Divergences: representation for the transverse set, arithmetic for the rest.

  dstat(1) = max(dstat(1), abs(real(x, rp) - x0(ip)))
  dstat(2) = max(dstat(2), abs(real(px, rp) - px0(ip)))
  dstat(3) = max(dstat(3), abs(real(y, rp) - y0(ip)))
  dstat(4) = max(dstat(4), abs(real(py, rp) - py0(ip)))
  dstat(5) = max(dstat(5), abs(real(pz, rp) - sl%pz(ip)))
  dstat(6) = max(dstat(6), abs(real(d_new, rp) - theta64))
enddo

! Scales: relative where a scale exists, absolute radians for theta.

sc(1) = maxval(abs(x0(1:n))) + 1e-30_rp
sc(2) = maxval(abs(px0(1:n))) + 1e-30_rp
sc(3) = maxval(abs(y0(1:n))) + 1e-30_rp
sc(4) = maxval(abs(py0(1:n))) + 1e-30_rp
sc(5) = maxval(abs(sl%pz(1:n))) + 1e-30_rp
sc(6) = 1

! The phasor scale is the slice's full charge weight, not |P64|: before bunching forms
! the phasor is a noise-level sum and a relative error against it means nothing. On
! this scale the row is the error in the bunching factor itself.

sc(7) = wsum + 1e-30_rp

fp32%div_slice(1:6, is) = dstat(1:6) / sc(1:6)
fp32%div_slice(7, is) = abs(cmplx(p32sum, kind=rp) - p64sum) / sc(7)

! The guard statistic: the median per-step phase-residual increment in ulps of the
! residual's own magnitude. The silent failure is this number reaching zero.

call median_inplace (incr, n, ulp_med)
fp32%ulp_slice(is) = ulp_med / max(real(spacing(maxval(abs(ks32 * t%dzr(1:n)))), rp), tiny(1.0_rp))

!------------------------------------------------------------------------------
contains

! The FP32 image of faw (fel_track_mod), the same expressions in sp.

function faw32 (xx, yy) result (value)
real(sp) xx, yy, value, ddx, ddy, ddt
ddx = xx - ax32
ddy = yy - ay32
if (st32 /= 0) then
  ddt = ct32 * ddx + st32 * ddy
  ddy = -st32 * ddx + ct32 * ddy
  ddx = ddt
endif
value = 1 + 0.5_sp * (kx32 * ddx*ddx + ky32 * ddy*ddy)
end function faw32

! The FP64 faw, for the comparison side (same expressions as fel_track_mod's).

function faw64 (xx, yy) result (value)
real(rp) xx, yy, value, ddx, ddy, ddt
ddx = xx - ax
ddy = yy - ay
if (sin_t /= 0) then
  ddt = cos_t * ddx + sin_t * ddy
  ddy = -sin_t * ddx + cos_t * ddy
  ddx = ddt
endif
value = 1 + 0.5_rp * (kx * ddx*ddx + ky * ddy*ddy)
end function faw64

! The FP32 gather: fel_grid_weights_pre's expressions in sp, the four corner values
! converted to FP32 before the blend, as a resident FP32 field would hold them.

subroutine gather32 (xx, yy, cp, ong)
real(sp) xx, yy, wwx, wwy
complex(sp) cp
logical ong
integer jx, jy
if (xx > -gmax32 .and. xx < gmax32 .and. yy > -gmax32 .and. yy < gmax32) then
  wwx = (xx + gmax32) / dxg32
  wwy = (yy + gmax32) / dyg32
  jx = int(floor(wwx));  jy = int(floor(wwy))
  wwx = 1 + real(jx, sp) - wwx
  wwy = 1 + real(jy, sp) - wwy
  jx = jx + 1;  jy = jy + 1

  ! FP32 rounding at the outer boundary can land one cell past where the FP64 form
  ! stops. That is an off-grid particle at this precision, and it counts as one.

  if (jx < 1 .or. jy < 1 .or. jx+1 > size(exfld,1) .or. jy+1 > size(exfld,2)) then
    cp = 0
    ong = .false.
    return
  endif
  cp =      cmplx(exfld(jx,   jy  ), kind=sp) * wwx * wwy
  cp = cp + cmplx(exfld(jx+1, jy  ), kind=sp) * (1-wwx) * wwy
  cp = cp + cmplx(exfld(jx,   jy+1), kind=sp) * wwx * (1-wwy)
  cp = cp + cmplx(exfld(jx+1, jy+1), kind=sp) * (1-wwx) * (1-wwy)
  ong = .true.
else
  cp = 0
  ong = .false.
endif
end subroutine gather32

! The FP32 image of fel_runge_kutta on (goff, delta), the same stage structure.

subroutine rk32 (g, d)
real(sp) g, d
real(sp) k2gg, k2pp, k3gg, k3pp, stpz
k2gg = 0;  k2pp = 0
call ode32 (g, d, k2gg, k2pp)
stpz = 0.5_sp * dz32
g = g + stpz * k2gg;  d = d + stpz * k2pp
k3gg = k2gg;  k3pp = k2pp
k2gg = 0;  k2pp = 0
call ode32 (g, d, k2gg, k2pp)
g = g + stpz * (k2gg - k3gg);  d = d + stpz * (k2pp - k3pp)
k3gg = k3gg / 6;  k3pp = k3pp / 6
k2gg = k2gg * (-0.5_sp);  k2pp = k2pp * (-0.5_sp)
call ode32 (g, d, k2gg, k2pp)
stpz = dz32
g = g + stpz * k2gg;  d = d + stpz * k2pp
k3gg = k3gg - k2gg;  k3pp = k3pp - k2pp
k2gg = k2gg * 2;  k2pp = k2pp * 2
call ode32 (g, d, k2gg, k2pp)
g = g + stpz * (k3gg + k2gg / 6.0_sp)
d = d + stpz * (k3pp + k2pp / 6.0_sp)
end subroutine rk32

! The FP32 image of fel_ode under the detuning reformulation (module header):
! k2pp gains 0.5*ks*(qres - q) in place of ks*(1 - 1/btpar0) + ku, whose sqrt
! argument is exactly 1 in FP32. b0 is 1 to FP32 resolution and is not formed.

subroutine ode32 (g, d, kg, kp)
real(sp) g, d, kg, kp
real(sp) s_t, c_t, gam_l, q_l, btper
complex(sp) ctmp
call sincos32 (d, s_t, c_t)
ctmp = rpart * (base32 * cmplx(c_t, -s_t, sp))
gam_l = gam032 + g
btper = btpar + (-2.0_sp / ks32) * real(ctmp, sp)
q_l = btper / (gam_l * gam_l)
kp = kp + 0.5_sp * ks32 * (qres32 - q_l)
kg = kg + aimag(ctmp) / gam_l
end subroutine ode32

! The FP32 pair, plain intrinsics on the small angle.

subroutine sincos32 (a, s, c)
real(sp) a, s, c
s = sin(a)
c = cos(a)
end subroutine sincos32

! The median of the first nn entries, scratch reordered in place. Heapsort: O(n log n),
! no recursion, no degenerate case. This runs per slice per step, so it must be cheap
! and it must terminate on any input. The sift is inlined twice because an internal
! procedure cannot contain one of its own.

subroutine median_inplace (v, nn, med)
real(rp) v(:), med, tmp
integer nn, i, iend, root, child

do i = nn/2, 1, -1
  root = i
  do while (2*root <= nn)
    child = 2*root
    if (child < nn) then
      if (v(child+1) > v(child)) child = child + 1
    endif
    if (v(root) >= v(child)) exit
    tmp = v(root);  v(root) = v(child);  v(child) = tmp
    root = child
  enddo
enddo

do iend = nn, 2, -1
  tmp = v(1);  v(1) = v(iend);  v(iend) = tmp
  root = 1
  do while (2*root <= iend - 1)
    child = 2*root
    if (child < iend - 1) then
      if (v(child+1) > v(child)) child = child + 1
    endif
    if (v(root) >= v(child)) exit
    tmp = v(root);  v(root) = v(child);  v(child) = tmp
    root = child
  enddo
enddo

if (mod(nn, 2) == 1) then
  med = v((nn+1)/2)
else
  med = 0.5_rp * (v(nn/2) + v(nn/2+1))
endif

end subroutine median_inplace

end subroutine fel_fp32_twin_slice

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_fp32_step_close (fp32, err_flag)
!
! Routine to reduce the step's per-slice divergences, write the row, and run the
! guard. Serial, after the caller's parallel loop. The guard refuses when the median
! per-step residual increment falls below 32 ulps of the residual: below that the
! FP32 representation is absorbing the physics, which is the silent failure the
! design brief records.
!-

subroutine fel_fp32_step_close (fp32, err_flag)

type (fel_fp32_struct) fp32
logical err_flag
real(rp) row(fel_fp32_nq$), gulp
integer iq
character(*), parameter :: r_name = 'fel_fp32_step_close'

!

err_flag = .false.
fp32%istep = fp32%istep + 1

do iq = 1, fel_fp32_nq$
  row(iq) = maxval(fp32%div_slice(iq, :))
  fp32%worst(iq) = max(fp32%worst(iq), row(iq))
enddo
gulp = minval(fp32%ulp_slice)
fp32%ulp_min = min(fp32%ulp_min, gulp)

write (fp32%iu, '(i8, 8es13.4)') fp32%istep, row, gulp

if (gulp < 32) then
  call out_io (s_error$, r_name, 'FP32 RESIDUAL GUARD: THE MEDIAN PER-STEP PHASE INCREMENT IS \es10.2\ ULPS', &
        'OF THE RESIDUAL, BELOW THE FLOOR OF 32. AT THIS GRANULARITY FP32 ABSORBS THE PHYSICS', &
        'SILENTLY AND A RUN LOOKS CLEAN WITH NO BUNCHING (doc/performance.md), SO IT REFUSES INSTEAD.', &
        r_array = [gulp])
  err_flag = .true.
endif

end subroutine fel_fp32_step_close

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_fp32_close (fp32)
!
! Routine to write the run summary block and close the stream.
!-

subroutine fel_fp32_close (fp32)

type (fel_fp32_struct) fp32
integer iq

!

if (.not. fp32%on) return
if (fp32%iu == 0) return

write (fp32%iu, '(a)') '#'
do iq = 1, fel_fp32_nq$
  write (fp32%iu, '(2a, es13.4)') 'worst_', trim(fel_fp32_qname(iq)), fp32%worst(iq)
enddo
write (fp32%iu, '(a, es13.4)') 'guard_ulp_min ', fp32%ulp_min
write (fp32%iu, '(a, i0)') 'steps ', fp32%istep
close (fp32%iu)
fp32%iu = 0

end subroutine fel_fp32_close

end module fel_fp32_mod
