!+
! Module fel_import_mod
!
! This module is the importdistribution equivalent: it resamples a Bmad bunch_struct --
! arbitrary times, arbitrary weights -- into the evenly spaced, equal-population slices
! the FEL step wants, by Genesis 1.3 v4's own method, transcribed from SDDSBeam.cpp
! (the class name is historical: it reads plain HDF5).
! Method, Genesis provenance, and validation: lucifer/doc/fel-physics.tex
! (sec:import). The only weighted generalization: the slice current is c*sum(w)/dslen
! where Genesis counts particles (identical for uniform weights).
!
! Longitudinal mapping: a bunch particle's window position is the port's own
! tau = -z/beta = c*(t - t_ref), min-shifted to zero exactly as Genesis min-shifts its
! file's s = -c*t, so both codes bin the identical particle set identically when the
! distribution file is written with t = -tau/c (sec:import).
!
! Genesis's match/center transforms are NOT ported, by decision: they exist because
! Genesis lattices carry no optics, so an imported bunch must be rematched by hand. A
! Bmad lattice carries its Twiss, and init_beam_distribution generates bunches matched
! to the lattice element already -- the transform would be a second way to say what
! the lattice says. (An openPMD bunch that genuinely needs rematching is a Bmad
! tracking problem upstream of the FEL, not an import option.) The bunch moments ARE
! still measured -- an UNWEIGHTED analysis in Genesis's analyse form, unweighted
! deliberately: coincident split-weight copies leave every moment bit-identical, which
! the invariance check relies on.
!
! Genesis quirks found by reading and NOT transcribed as functional (sec:import): the
! align/align_start/align_end parameters are parsed but never used in v4; the
! shotnoise flag is read but never consulted -- the import applies noise
! unconditionally, skipping only zero-current slices. Both behaviors are kept as
! Genesis has them (noise always on for nonzero current). one4one is out of scope:
! per-particle weights supersede it.
!
! The RNG is Bmad's (ran_uniform), NOT a transcription of Genesis's RandomU, so
! everything the RNG touches is validated statistically; the current profile and the
! analyse/match/center moments are RNG-free and check exactly.
!-

module fel_import_mod

use fel_beam_mod

implicit none

!+
! Struct fel_import_param_struct
!
! The &importdistribution knobs, named after Genesis's where one exists.
!-

type fel_import_param_struct
  real(rp) :: slicewidth = 0.01_rp       ! Sampling window / bunch length (Genesis ds).
  integer :: npart = 8192                ! Macroparticles per slice after resampling.
  integer :: nbins = 4                   ! Beamlet size (quiet-load bins).
  integer :: nslice = 0                  ! 0: round(bunch_length/slice_spacing), Genesis's rule.
end type

private analyse_window

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_import_bunch (bunch, gamma0, lambda0, slice_spacing, prm, fbeam, err_flag, moments_out)
!
! Routine to resample a bunch_struct into FEL slices by the transcribed Genesis method
! (module header). gamma0 is the FEL reference, resolved from the lattice by the caller
! (it is also the filler energy for empty slices); lambda0 the radiation wavelength;
! slice_spacing = sample*lambda0. The caller seeds Bmad's RNG.
!
! Input:
!   bunch           -- bunch_struct: Bunch to resample; dead particles are skipped.
!   gamma0          -- real(rp): FEL reference gamma; sets fbeam%p0c.
!   lambda0         -- real(rp): Radiation wavelength [m].
!   slice_spacing   -- real(rp): Longitudinal slice spacing, sample*lambda0 [m].
!   prm             -- fel_import_param_struct: The &importdistribution knobs.
!
! Output:
!   fbeam           -- fel_beam_struct: The resampled packed beam.
!   err_flag        -- logical: Set True on error, False otherwise.
!   moments_out(11) -- real(rp), optional: The RNG-free analysis moments in Genesis's
!                       analyse order (gavg, xavg, pxavg, yavg, pyavg, ex, ey, bx, by,
!                       ax, ay) -- the deterministic quantities the exactness checks read.
!-

subroutine fel_import_bunch (bunch, gamma0, lambda0, slice_spacing, prm, fbeam, err_flag, moments_out)

type (bunch_struct), target :: bunch
type (fel_import_param_struct) prm
type (fel_beam_struct), target :: fbeam
type (fel_slice_struct), pointer :: sl
type (coord_struct), pointer :: cp

real(rp) gamma0, lambda0, slice_spacing
real(rp), optional :: moments_out(11)

! The whole-bunch working set, in Genesis's variables: s (window position), gamma,
! x, y and the SLOPES xp, yp. Weight rides along for the current sums.
real(rp), allocatable :: s(:), gam(:), x(:), y(:), xp(:), yp(:), wt(:)
! The per-slice candidate set, in Genesis's internal Particle coordinates
! (px = xp*gamma after the slope conversion).
real(rp), allocatable :: cg(:), cx(:), cy(:), cpx(:), cpy(:), theta(:), wk(:)

real(rp) mom(11), gavg, xavg, pxavg, yavg, pyavg, ex, ey, bx, by, ax, ay
real(rp) smin, ttotal, dslen, sloc, p_mc, p0_mc, u, ks_l, cur, ne, w_part, beta
integer n, nalive, i, ip, islice, nslice, mpart, ncand, n_clamp
logical err_flag

character(*), parameter :: r_name = 'fel_import_bunch'

!

err_flag = .true.

if (prm%npart < 1 .or. prm%nbins < 1 .or. mod(prm%npart, prm%nbins) /= 0) then
  call out_io (s_error$, r_name, 'NPART MUST BE A POSITIVE MULTIPLE OF NBINS (Genesis''s rule).')
  return
endif
if (prm%slicewidth <= 0) then
  call out_io (s_error$, r_name, 'SLICEWIDTH MUST BE POSITIVE.')
  return
endif

! The bunch in Genesis's variables. tau = -z/beta = c*(t - t_ref), the same chart theta
! derives from everywhere else in this port; slopes xp = Px/P = vec(2)/(1+pz), an exact
! round trip of the Genesis reconstruction px = xp*gamma (gamma*beta_x approx) through
! this port's stored px = Px/p0.

n = size(bunch%particle)
nalive = count(bunch%particle%state == alive$)
if (nalive < 1) then
  call out_io (s_error$, r_name, 'BUNCH HAS NO LIVE PARTICLES.')
  return
endif

allocate (s(nalive), gam(nalive), x(nalive), y(nalive), xp(nalive), yp(nalive), wt(nalive))

i = 0
do ip = 1, n
  cp => bunch%particle(ip)
  if (cp%state /= alive$) cycle
  i = i + 1
  p_mc = (1 + cp%vec(6)) * cp%p0c / m_electron
  gam(i) = sqrt(p_mc**2 + 1)
  s(i) = -cp%vec(5) * gam(i) / p_mc               ! tau = -z/beta
  x(i) = cp%vec(1);  y(i) = cp%vec(3)
  xp(i) = cp%vec(2) / (1 + cp%vec(6))
  yp(i) = cp%vec(4) / (1 + cp%vec(6))
  wt(i) = cp%charge
enddo

! A bunch with no charge imports as a perfectly dark beam -- every window current
! zero, every weight zero, a run that tracks and produces nothing, silently. The
! usual causes are an openPMD file without charge data and an unset
! beam_init%bunch_charge. Refuse by name.

if (sum(wt) <= 0) then
  call out_io (s_error$, r_name, 'BUNCH HAS ZERO TOTAL CHARGE; NOTHING WOULD LASE.', &
    'AN openPMD FILE WITHOUT CHARGE DATA, OR AN UNSET beam_init%bunch_charge, IMPORTS DARK.')
  return
endif

smin = minval(s);  s = s - smin                    ! Genesis's min shift (sec:import).
ttotal = maxval(s)
if (ttotal <= 0) then
  call out_io (s_error$, r_name, 'BUNCH HAS ZERO LENGTH; NOTHING TO SLICE.')
  return
endif

nslice = prm%nslice
if (nslice < 1) nslice = max(1, nint(ttotal / slice_spacing))   ! Genesis's rule (sec:import).

! The bunch moments, in Genesis's analyse form (unweighted, strict window bounds --
! the two extreme particles are excluded, exactly as Genesis's eval defaults do).

call analyse_window (s, gam, x, y, xp, yp, 0.0_rp, ttotal, &
                     gavg, xavg, pxavg, yavg, pyavg, ex, ey, bx, by, ax, ay)

mom = [gavg, xavg, pxavg, yavg, pyavg, ex, ey, bx, by, ax, ay]
if (present(moments_out)) moments_out = mom

! The beam container.

p0_mc = sqrt(gamma0**2 - 1)
fbeam%p0c = p0_mc * m_electron
fbeam%phi0 = 0
fbeam%wavelength = lambda0
fbeam%slice_spacing = slice_spacing
fbeam%s0 = 0
fbeam%nbins = prm%nbins
fbeam%one4one = .false.
if (allocated(fbeam%slice)) deallocate(fbeam%slice)
allocate (fbeam%slice(nslice))

! The slice loop (fel-physics.tex sec:import). Slice centers at
! (islice-1)*slice_spacing; candidates from the dslen window, strict inequalities.

dslen = prm%slicewidth * ttotal
mpart = prm%npart / prm%nbins
ks_l = twopi / lambda0
n_clamp = 0

allocate (cg(prm%npart), cx(prm%npart), cy(prm%npart), cpx(prm%npart), cpy(prm%npart), &
          theta(prm%npart), wk(prm%npart))

do islice = 1, nslice
  sloc = (islice - 1) * slice_spacing

  ! Candidates: window membership and the weighted current from the same window.
  ncand = 0
  cur = 0
  do i = 1, nalive
    if (s(i) > sloc - 0.5_rp*dslen .and. s(i) < sloc + 0.5_rp*dslen) then
      cur = cur + wt(i)
      if (ncand < mpart) then          ! Keep up to mpart directly; excess handled below.
        ncand = ncand + 1
        cg(ncand) = gam(i); cx(ncand) = x(i); cy(ncand) = y(i)
        cpx(ncand) = xp(i) * gam(i)    ! Genesis's slope-to-momentum conversion.
        cpy(ncand) = yp(i) * gam(i)
      else
        ! Genesis copies every candidate then deletes at random (removeParticles:
        ! overwrite a random index with the last, shrink). Deleting-by-random-index
        ! from the full set is equivalent to reservoir-style replacement only if the
        ! transcription keeps Genesis's exact rule, so transcribe the rule: append,
        ! then delete. Appending needs storage: do it in two passes instead.
        ncand = ncand + 1              ! Count only; second pass below when oversized.
      endif
    endif
  enddo
  cur = cur * c_light / dslen                       ! c*sum(w)/dslen

  if (ncand > mpart) then
    ! Second pass, transcribing removeParticles faithfully: gather ALL candidates,
    ! then repeatedly overwrite a random index with the last and shrink.
    call gather_and_remove ()
  endif

  call fill_slice ()
enddo

if (n_clamp > 0) then
  call out_io (s_warn$, r_name, 'SOME BEAMLET NOISE DRAWS HAD FEWER THAN ONE REAL ELECTRON', &
    '(CLAMPED TO 1, AS GENESIS DOES SILENTLY); THE NOISE THERE IS NOT PHYSICAL.')
endif

err_flag = .false.

!------------------------------------------------------------------------------
contains

!+
! Subroutine gather_and_remove ()
!
! Routine to gather every window candidate (needs its own storage: there can be far more
! than npart), then apply Genesis's removeParticles: overwrite a random index with the
! last and shrink until mpart remain.
!-

subroutine gather_and_remove ()

real(rp), allocatable :: ag(:), axx(:), ayy(:), apx(:), apy(:)
integer m, k
real(rp) uu

allocate (ag(ncand), axx(ncand), ayy(ncand), apx(ncand), apy(ncand))
m = 0
do k = 1, nalive
  if (s(k) > sloc - 0.5_rp*dslen .and. s(k) < sloc + 0.5_rp*dslen) then
    m = m + 1
    ag(m) = gam(k); axx(m) = x(k); ayy(m) = y(k)
    apx(m) = xp(k) * gam(k); apy(m) = yp(k) * gam(k)
  endif
enddo

do while (m > mpart)                               ! removeParticles's exact rule.
  call ran_uniform (uu)
  k = int(m * uu) + 1
  if (k > m) k = m
  ag(k) = ag(m); axx(k) = axx(m); ayy(k) = ayy(m)
  apx(k) = apx(m); apy(k) = apy(m)
  m = m - 1
enddo

cg(1:mpart) = ag(1:mpart); cx(1:mpart) = axx(1:mpart); cy(1:mpart) = ayy(1:mpart)
cpx(1:mpart) = apx(1:mpart); cpy(1:mpart) = apy(1:mpart)
ncand = mpart

end subroutine gather_and_remove

!..............................................................................
!+
! Subroutine fill_slice ()
!
! Routine to bring the candidate set to mpart seeds (Genesis's addParticles when short),
! refill theta, mirror into nbins beamlet copies, impose the shot noise, and store the
! slice.
!-

subroutine fill_slice ()

real(rp) g1, x1, y1, px1, py1, g2, x2, y2, px2, py2, scl, rmin, r, tmp, uu
real(rp) gam_p, beta_p
integer nd, nd0, k, j, n1, n2, i1, i2

! addParticles (sec:import): empty -> one filler at the reference energy;
! singleton -> mirror it; then interpolate up to mpart.

nd = ncand
if (nd == 0) then
  cg(1) = gamma0; cx(1) = 0; cy(1) = 0; cpx(1) = 0; cpy(1) = 0
  nd = 1
endif
if (nd == 1 .and. mpart > 1) then
  cg(2) = cg(1); cx(2) = cx(1); cy(2) = cy(1); cpx(2) = cpx(1); cpy(2) = cpy(1)
  nd = 2
endif

if (nd < mpart) then
  ! Normalize the five coordinates to zero mean, unit rms (cold dimensions keep
  ! scale 1, as Genesis's).
  g1 = sum(cg(1:nd))/nd;   g2 = sqrt(abs(sum(cg(1:nd)**2)/nd - g1**2))
  x1 = sum(cx(1:nd))/nd;   x2 = sqrt(abs(sum(cx(1:nd)**2)/nd - x1**2))
  y1 = sum(cy(1:nd))/nd;   y2 = sqrt(abs(sum(cy(1:nd)**2)/nd - y1**2))
  px1 = sum(cpx(1:nd))/nd; px2 = sqrt(abs(sum(cpx(1:nd)**2)/nd - px1**2))
  py1 = sum(cpy(1:nd))/nd; py2 = sqrt(abs(sum(cpy(1:nd)**2)/nd - py1**2))
  if (g2 == 0)  then; g2 = 1;  else; g2 = 1/g2;   endif
  if (x2 == 0)  then; x2 = 1;  else; x2 = 1/x2;   endif
  if (y2 == 0)  then; y2 = 1;  else; y2 = 1/y2;   endif
  if (px2 == 0) then; px2 = 1; else; px2 = 1/px2; endif
  if (py2 == 0) then; py2 = 1; else; py2 = 1/py2; endif
  cg(1:nd) = (cg(1:nd) - g1) * g2
  cx(1:nd) = (cx(1:nd) - x1) * x2
  cy(1:nd) = (cy(1:nd) - y1) * y2
  cpx(1:nd) = (cpx(1:nd) - px1) * px2
  cpy(1:nd) = (cpy(1:nd) - py1) * py2

  nd0 = nd                                         ! Parents only from the originals.
  do while (nd < mpart)
    call ran_uniform (uu)
    n1 = int(nd0 * uu) + 1
    if (n1 > nd0) n1 = nd0
    rmin = 1e9_rp
    n2 = n1
    do k = 1, nd0                                  ! Nearest neighbor, random metric.
      if (k == n1) cycle                           ! Genesis's distance skips self.
      tmp = cg(n1)-cg(k);  call ran_uniform (uu);  r = tmp*tmp*uu
      tmp = cx(n1)-cx(k);  call ran_uniform (uu);  r = r + tmp*tmp*uu
      tmp = cy(n1)-cy(k);  call ran_uniform (uu);  r = r + tmp*tmp*uu
      tmp = cpx(n1)-cpx(k); call ran_uniform (uu); r = r + tmp*tmp*uu
      tmp = cpy(n1)-cpy(k); call ran_uniform (uu); r = r + tmp*tmp*uu
      if (r < rmin) then
        n2 = k
        rmin = r
      endif
    enddo
    nd = nd + 1
    call ran_uniform (uu)
    cg(nd) = 0.5_rp*(cg(n1)+cg(n2)) + (2*uu-1)*(cg(n1)-cg(n2))
    call ran_uniform (uu)
    cx(nd) = 0.5_rp*(cx(n1)+cx(n2)) + (2*uu-1)*(cx(n1)-cx(n2))
    call ran_uniform (uu)
    cpx(nd) = 0.5_rp*(cpx(n1)+cpx(n2)) + (2*uu-1)*(cpx(n1)-cpx(n2))
    call ran_uniform (uu)
    cy(nd) = 0.5_rp*(cy(n1)+cy(n2)) + (2*uu-1)*(cy(n1)-cy(n2))
    call ran_uniform (uu)
    cpy(nd) = 0.5_rp*(cpy(n1)+cpy(n2)) + (2*uu-1)*(cpy(n1)-cpy(n2))
  enddo

  cg(1:nd) = cg(1:nd)/g2 + g1                      ! Scale back.
  cx(1:nd) = cx(1:nd)/x2 + x1
  cy(1:nd) = cy(1:nd)/y2 + y1
  cpx(1:nd) = cpx(1:nd)/px2 + px1
  cpy(1:nd) = cpy(1:nd)/py2 + py1
endif

! theta refilled completely new over one beamlet spacing (sec:import), then the
! beamlet mirroring: seed i lands at nbins consecutive indices.

do k = 1, mpart
  call ran_uniform (uu)
  wk(k) = (twopi / prm%nbins) * uu
enddo
do k = mpart, 1, -1
  i1 = k
  i2 = prm%nbins * (k-1)
  do j = 0, prm%nbins - 1
    cg(i2+j+1) = cg(i1); cx(i2+j+1) = cx(i1); cy(i2+j+1) = cy(i1)
    cpx(i2+j+1) = cpx(i1); cpy(i2+j+1) = cpy(i1)
    theta(i2+j+1) = wk(i1) + j * twopi / prm%nbins
  enddo
enddo

! Uniform per-slice weights from the window current (Genesis's dQ semantics), the
! shared Fawley noise with ne = charge/e (skipped for empty slices, as Genesis's),
! then the stored chart: z = beta*theta/ks, px = (gamma*beta_x)/p0_mc.

w_part = cur * slice_spacing / (c_light * prm%npart)
ne = anint(cur * slice_spacing / (e_charge * c_light))

sl => fbeam%slice(islice)
call fel_slice_reallocate (sl, prm%npart)
sl%n = prm%npart
sl%weight(1:prm%npart) = w_part

if (ne > 0) then
  call fel_fawley_noise (theta(1:prm%npart), sl%weight(1:prm%npart), prm%npart, prm%nbins, n_clamp)
endif

do k = 1, prm%npart
  gam_p = cg(k)
  beta_p = sqrt(gam_p**2 - 1) / gam_p
  sl%x(k) = cx(k);  sl%px(k) = cpx(k) / p0_mc
  sl%y(k) = cy(k);  sl%py(k) = cpy(k) / p0_mc
  sl%z(k) = beta_p * theta(k) / ks_l
  sl%pz(k) = (sqrt(gam_p**2 - 1) - p0_mc) / p0_mc
enddo

end subroutine fill_slice

end subroutine fel_import_bunch

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine analyse_window (s, gam, x, y, xp, yp, s0, s1, gavg, xavg, pxavg, yavg, pyavg, ex, ey, bx, by, ax, ay)
!
! Routine to compute Genesis's analyse moments (fel-physics.tex sec:import): UNWEIGHTED
! means, variances and Twiss over the particles with s0 < s < s1, on the slopes.
! Emittance ex = sqrt(|var_x*var_px - cov^2|)*gavg (a normalized emittance through the
! mean energy), bx = var_x*gavg/ex, ax = -cov*gavg/ex. Unweighted deliberately -- see
! the module header.
!
! Input:
!   s(:)          -- real(rp): Window positions tau [m], min-shifted to zero.
!   gam(:)        -- real(rp): Lorentz factors.
!   x(:), y(:)    -- real(rp): Transverse positions [m].
!   xp(:), yp(:)  -- real(rp): Transverse slopes.
!   s0, s1        -- real(rp): Window bounds; only particles with s0 < s < s1 count
!                     (strict, so the two extreme particles are excluded, as Genesis's
!                     eval defaults do).
!
! Output:
!   gavg          -- real(rp): Mean gamma.
!   xavg, yavg    -- real(rp): Mean positions [m].
!   pxavg, pyavg  -- real(rp): Mean slopes.
!   ex, ey        -- real(rp): Normalized emittances through the mean energy.
!   bx, by        -- real(rp): Twiss beta functions [m].
!   ax, ay        -- real(rp): Twiss alpha functions.
!-

subroutine analyse_window (s, gam, x, y, xp, yp, s0, s1, gavg, xavg, pxavg, yavg, pyavg, ex, ey, bx, by, ax, ay)

real(rp) s(:), gam(:), x(:), y(:), xp(:), yp(:)
real(rp) s0, s1, gavg, xavg, pxavg, yavg, pyavg, ex, ey, bx, by, ax, ay
real(rp) a1, a2, b1, b2, ab, c1, c2, d1, d2, cd, e1, scl
integer i, ncount

!

ncount = 0
a1 = 0; a2 = 0; b1 = 0; b2 = 0; ab = 0
c1 = 0; c2 = 0; d1 = 0; d2 = 0; cd = 0
e1 = 0

do i = 1, size(s)
  if (s(i) <= s0 .or. s(i) >= s1) cycle
  ncount = ncount + 1
  a1 = a1 + x(i);   a2 = a2 + x(i)**2
  b1 = b1 + xp(i);  b2 = b2 + xp(i)**2
  ab = ab + x(i)*xp(i)
  c1 = c1 + y(i);   c2 = c2 + y(i)**2
  d1 = d1 + yp(i);  d2 = d2 + yp(i)**2
  cd = cd + y(i)*yp(i)
  e1 = e1 + gam(i)
enddo

if (ncount > 0) then
  scl = 1.0_rp / ncount
  gavg = e1*scl
  xavg = a1*scl;  pxavg = b1*scl
  yavg = c1*scl;  pyavg = d1*scl
  a2 = a2*scl; b2 = b2*scl; ab = ab*scl
  c2 = c2*scl; d2 = d2*scl; cd = cd*scl
else
  gavg = 0; xavg = 0; pxavg = 0; yavg = 0; pyavg = 0
  a2 = 0; b2 = 0; ab = 0; c2 = 0; d2 = 0; cd = 0
endif

ex = sqrt(abs((a2 - xavg**2)*(b2 - pxavg**2) - (ab - xavg*pxavg)**2)) * gavg
ey = sqrt(abs((c2 - yavg**2)*(d2 - pyavg**2) - (cd - yavg*pyavg)**2)) * gavg
bx = (a2 - xavg**2) / ex * gavg
by = (c2 - yavg**2) / ey * gavg
ax = -(ab - xavg*pxavg) * gavg / ex
ay = -(cd - yavg*pyavg) * gavg / ey

end subroutine analyse_window

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_write_genesis4_distribution (bunch, file_name, err_flag)
!
! Routine to write a bunch_struct as a Genesis 1.3 v4 DISTRIBUTION file (the
! &importdistribution input; not a dump): flat datasets t [s], p [gamma*beta], x, y [m],
! xp, yp [slopes], plus the total charge. t = -tau/c with tau = -z/beta, so Genesis's
! s = -c*t reproduces this port's window position exactly and both codes bin the
! identical particle set identically (fel-physics.tex sec:import). Dead particles are
! skipped.
!
! Input:
!   bunch       -- bunch_struct: Bunch to write.
!   file_name   -- character(*): File to create.
!
! Output:
!   err_flag    -- logical: Set True on error, False otherwise.
!-

subroutine fel_write_genesis4_distribution (bunch, file_name, err_flag)

type (bunch_struct), target :: bunch
type (coord_struct), pointer :: cp

integer(hid_t) f_id
integer ip, i, nalive, h5_err
real(rp) p_mc, g
real(rp), allocatable :: t(:), p(:), x(:), y(:), xp(:), yp(:)
logical err_flag, err
character(*) file_name
character(*), parameter :: r_name = 'fel_write_genesis4_distribution'

!

err_flag = .true.

nalive = count(bunch%particle%state == alive$)
if (nalive < 1) then
  call out_io (s_error$, r_name, 'BUNCH HAS NO LIVE PARTICLES.')
  return
endif

allocate (t(nalive), p(nalive), x(nalive), y(nalive), xp(nalive), yp(nalive))

i = 0
do ip = 1, size(bunch%particle)
  cp => bunch%particle(ip)
  if (cp%state /= alive$) cycle
  i = i + 1
  p_mc = (1 + cp%vec(6)) * cp%p0c / m_electron
  g = sqrt(p_mc**2 + 1)
  t(i) = cp%vec(5) * g / (p_mc * c_light)          ! -tau/c, tau = -z/beta
  p(i) = p_mc
  x(i) = cp%vec(1);  y(i) = cp%vec(3)
  xp(i) = cp%vec(2) / (1 + cp%vec(6))
  yp(i) = cp%vec(4) / (1 + cp%vec(6))
enddo

call hdf5_open_file (file_name, 'WRITE', f_id, err);  if (err) return

call hdf5_write_dataset_real (f_id, 'charge', [sum(bunch%particle%charge, mask = bunch%particle%state == alive$)], err)
if (err) return
call hdf5_write_dataset_real (f_id, 't',  t,  err);  if (err) return
call hdf5_write_dataset_real (f_id, 'p',  p,  err);  if (err) return
call hdf5_write_dataset_real (f_id, 'x',  x,  err);  if (err) return
call hdf5_write_dataset_real (f_id, 'y',  y,  err);  if (err) return
call hdf5_write_dataset_real (f_id, 'xp', xp, err);  if (err) return
call hdf5_write_dataset_real (f_id, 'yp', yp, err);  if (err) return

call H5Fclose_f (f_id, h5_err)

err_flag = .false.

end subroutine fel_write_genesis4_distribution

end module fel_import_mod
