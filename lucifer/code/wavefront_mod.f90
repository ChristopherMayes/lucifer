!+
! Module wavefront_mod
!
! Representation and free-space propagation of a paraxial radiation wavefront.
!
! A wavefront_struct holds the complex transverse electric field components Ex and Ey
! sampled on a uniform 3D grid of shape (nx, ny, nz) in V/m. The first two indices are
! the transverse coordinates, the third is the longitudinal slice index. The wavefront
! propagates in the +z direction with the pulse head at max(z).
!
! This mirrors the Wavefront class of openPMD-beamphysics
! (beamphysics/wavefront/wavefront.py) closely enough that the two can be compared field
! by field. See lucifer/wavefront/tests/run_validation.sh, which drifts the same input
! through both and reports the largest relative difference.
!
! openPMD EXT_Wavefront input and output for this structure is in
! wavefront_openpmd_mod.
!
! Index order. Fortran is column major, so Ex(:,:,iz) (one longitudinal slice) is
! contiguous. That is what a per-slice 2D FFTW plan wants, and what parallelising over
! slices wants. The Python class uses the same (nx, ny, nz) shape but in C order, so
! there its slices are strided.
!
! Grid convention, matching the Python class: the transverse grids are centered on zero,
! with x running from -(nx-1)*dx/2 to +(nx-1)*dx/2, and likewise for y. The slice axis
! runs from the tail at iz = 1 to the head at iz = nz. A whole window's first slice sits
! at z = 0. There is no separate origin offset. The one absolute position carried is
! s_position, the beamline position of the plane, which is openPMD's zCoordinate.
! Genesis's refposition is a different quantity, the in-bunch coordinate of a file's
! first slice. This structure does not carry it. A writer of that format computes it
! from the window (fel_io_mod).
!-

module wavefront_mod

use, intrinsic :: iso_c_binding
use sim_utils

use, intrinsic :: ieee_arithmetic

implicit none

! Kind of the complex field arrays.
!
! This is the one point of control for the field precision. Every field array, and every
! routine argument carrying field values, is declared complex(wf_rp). Changing this
! parameter changes the representation everywhere without touching a call site.
!
! Note: wf_rp is rp (double). Changing it alone does not give a single precision path.
! The plan cache below binds the double precision FFTW entry points, and its work buffer
! is complex(wf_rp). So a single precision field fails to compile at fftw_execute_dft
! with a type mismatch. That is intended. Transforming a single precision field in double
! precision would be a silent fallback. The compile error says that the fftwf_* entry
! points still have to be bound.

integer, parameter :: wf_rp = rp

! Transform directions for wavefront_fft2. The sign convention is that of FFTW and
! numpy.fft. Forward carries exp(-i k x) and backward carries exp(+i k x). Neither
! direction applies the 1/(nx*ny) factor. The caller does.
!
! Note: Interchanging these two values does not change what wavefront_drift computes, and
! no test can detect it. The backward transform is the forward transform with k negated.
! The propagation kernel depends on kx^2 + ky^2 and so is even in k. Substituting k -> -k
! in the sum over k recovers the original expression term for term. So the composition
! forward, multiply by the kernel, backward is invariant under the interchange. The
! individual signs matter for any operator whose kernel is not even in k.

integer, parameter :: wf_fft_forward$ = -1, wf_fft_backward$ = 1

!+
! Struct wavefront_struct
!
! Ex, Ey are allocated independently. An unallocated component means that polarisation
! component is absent. That is how the Python class's Ex = None is represented. At least
! one of the two must be allocated for the wavefront to be usable. wavefront_check tests
! that, along with the rest of the Python class's __post_init__ validation.
!-

type wavefront_struct
  complex(wf_rp), allocatable :: Ex(:,:,:)   ! x polarized field [V/m], shape (nx, ny, nz).
  complex(wf_rp), allocatable :: Ey(:,:,:)   ! y polarized field [V/m], shape (nx, ny, nz).
  real(rp) :: dx = 1                         ! Transverse grid spacing in x [m].
  real(rp) :: dy = 1                         ! Transverse grid spacing in y [m].
  real(rp) :: dz = 1                         ! Slice spacing [m].
  real(rp) :: wavelength = 1                 ! Central wavelength [m].
  real(rp) :: s_position = 0                 ! Beamline position of the plane [m]:
                                             !   openPMD's zCoordinate. wavefront_drift
                                             !   advances it by the drift, and the FEL walk
                                             !   sets it to the lattice position it holds.
                                             !   Never Genesis's refposition, which is a
                                             !   coordinate within the pulse.
end type

!+
! Struct wavefront_params_struct
!
! Summary statistics of one field slice. The field analog of Bmad's bunch_params_struct,
! and named to match it. centroid and sigma are the intensity-weighted first and second
! Wigner moments over the transverse phase space (x, theta_x, y, theta_y). So sizes are
! sqrt(sigma(1,1)) in both structs, and free-space propagation is the ABCD map on sigma.
! sigma_x^2(z) is quadratic in z. That is how banked slices are folded into pulse
! statistics without numerical propagation. emit = sqrt(det of a plane's 2x2 block) is
! the field-quality analog (= M^2 lambda/4pi). Pulse-level values are pooled from slice
! instances downstream (energy-weighted mean of sigmas plus variance of centroids) and
! are never stored. Units are fixed Bmad units.
!
! The theta rows of sigma cost FFT-side sums, so they are filled only where needed
! (element ends, bank time). angle_moments_valid says whether they were, as twiss_valid
! does.
!-

type wavefront_params_struct
  real(rp) :: centroid(4) = 0        ! Intensity-weighted (x, theta_x, y, theta_y) [m, rad].
  real(rp) :: sigma(4,4) = 0         ! Second moments about the centroid [m^2, m rad, rad^2].
  real(rp) :: energy = 0             ! Field energy of the slice [J].
  real(rp) :: power = 0              ! Slice power [W].
  real(rp) :: on_axis_intensity = 0  ! Intensity at the grid center [W/m^2].
  real(rp) :: emit_x = 0, emit_y = 0 ! sqrt(det sigma_plane) [m rad].
  real(rp) :: s = -1                 ! Longitudinal position of evaluation [m].
  logical :: angle_moments_valid = .false.
end type

! Cached FFTW plan state for wavefront_fft2, private to the module, one cache per thread.
!
! This follows bmad/space_charge/fft_interface_mod.f90: plans are created once per
! transverse grid size and reused. There are two differences. First, the work buffer is
! allocated by fftw_alloc_complex, and transforms run on it rather than on the caller's
! array. So the plan always sees the alignment it was created with. FFTW's new-array
! execute rule makes executing a plan on a differently aligned array undefined.
! Second, the cache is threadprivate. Every OpenMP thread carries its own plans and work
! buffer, so wavefront_fft2 can be called from a parallel loop over slices.
! wavefront_drift and the FEL field solve call it that way. The cost is one planner run
! and one nx*ny buffer per thread, once per grid size.
!
! The critical section inside wavefront_fft2_plan guards FFTW's rule that plan creation
! is not reentrant, and nothing else. The FFTW planner is globally serialised however
! many per-thread buffers exist. Plan execution is thread safe by FFTW's own guarantee,
! and each execution here touches only the calling thread's buffer.

! The transform itself is single threaded. The parallel axis for a wavefront is the slice
! index, not the transverse transform. So fftw_plan_with_nthreads is not called.
! Threading a small 2D transform under a parallel loop over slices would oversubscribe.

type(C_PTR), private, save :: wf_plan_fwd = C_NULL_PTR
type(C_PTR), private, save :: wf_plan_bwd = C_NULL_PTR
type(C_PTR), private, save :: wf_buf_cptr = C_NULL_PTR
complex(wf_rp), private, pointer, save :: wf_buf(:,:) => null()
integer, private, save :: wf_cache_nx = 0, wf_cache_ny = 0
!$OMP threadprivate(wf_plan_fwd, wf_plan_bwd, wf_buf_cptr, wf_buf, wf_cache_nx, wf_cache_ny)

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_init (wf, nx, ny, nz, dx, dy, dz, wavelength, polarization, s_position)
!
! Routine to allocate the field arrays of a wavefront and set its grid parameters.
! The allocated components are zeroed.
!
! Input:
!   nx, ny, nz    -- integer: Grid dimensions.
!   dx, dy, dz    -- real(rp): Grid spacings [m]. Must be positive.
!   wavelength    -- real(rp): Central wavelength [m]. Must be positive.
!   polarization  -- character(*), optional: Which components to allocate. One of 'x',
!                      'y' or 'xy'. Default is 'x'.
!   s_position    -- real(rp), optional: Beamline position of the plane [m]. Default is 0.
!
! Output:
!   wf            -- wavefront_struct: Initialized wavefront.
!-

subroutine wavefront_init (wf, nx, ny, nz, dx, dy, dz, wavelength, polarization, s_position)

type (wavefront_struct) wf
integer nx, ny, nz
real(rp) dx, dy, dz, wavelength
real(rp), optional :: s_position
character(*), optional :: polarization
character(*), parameter :: r_name = 'wavefront_init'
character(4) pol

!

pol = string_option('x', polarization)

if (pol /= 'x' .and. pol /= 'y' .and. pol /= 'xy') then
  call out_io (s_error$, r_name, 'POLARIZATION MUST BE "x", "y" OR "xy". GOT: ' // pol)
  return
endif

if (allocated(wf%Ex)) deallocate(wf%Ex)
if (allocated(wf%Ey)) deallocate(wf%Ey)

if (index(pol, 'x') /= 0) then
  allocate (wf%Ex(nx, ny, nz))
  wf%Ex = 0
endif

if (index(pol, 'y') /= 0) then
  allocate (wf%Ey(nx, ny, nz))
  wf%Ey = 0
endif

wf%dx = dx
wf%dy = dy
wf%dz = dz
wf%wavelength = wavelength
wf%s_position = real_option(0.0_rp, s_position)

end subroutine wavefront_init

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_check (wf, err_flag)
!
! Routine to validate a wavefront, the Fortran counterpart of the Python class's
! __post_init__. At least one field component must be present. The two components must
! have the same shape if both are present. dx, dy, dz and wavelength must be positive.
!
! Input:
!   wf        -- wavefront_struct: Wavefront to check.
!
! Output:
!   err_flag  -- logical: Set True if the wavefront is not valid, False otherwise.
!-

subroutine wavefront_check (wf, err_flag)

type (wavefront_struct) wf
logical err_flag
character(*), parameter :: r_name = 'wavefront_check'

!

err_flag = .true.

if (.not. allocated(wf%Ex) .and. .not. allocated(wf%Ey)) then
  call out_io (s_error$, r_name, 'AT LEAST ONE OF Ex OR Ey MUST BE PRESENT IN A WAVEFRONT.')
  return
endif

if (allocated(wf%Ex) .and. allocated(wf%Ey)) then
  if (any(shape(wf%Ex) /= shape(wf%Ey))) then
    call out_io (s_error$, r_name, 'Ex SHAPE \3i6\ DOES NOT MATCH Ey SHAPE \3i6\.', &
                                                 i_array = [shape(wf%Ex), shape(wf%Ey)])
    return
  endif
endif

if (wf%dx <= 0 .or. wf%dy <= 0 .or. wf%dz <= 0 .or. wf%wavelength <= 0) then
  call out_io (s_error$, r_name, 'dx, dy, dz AND wavelength MUST ALL BE POSITIVE. GOT: \4es14.6\', &
                                                 r_array = [wf%dx, wf%dy, wf%dz, wf%wavelength])
  return
endif

err_flag = .false.

end subroutine wavefront_check

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_shape (wf) result (n_grid)
!
! Routine to return the shape of the field arrays of a wavefront.
!
! Input:
!   wf          -- wavefront_struct: Wavefront.
!
! Output:
!   n_grid(3)   -- integer: [nx, ny, nz]. Set to zero if neither component is present.
!-

function wavefront_shape (wf) result (n_grid)

type (wavefront_struct) wf
integer n_grid(3)

!

if (allocated(wf%Ex)) then
  n_grid = shape(wf%Ex)
elseif (allocated(wf%Ey)) then
  n_grid = shape(wf%Ey)
else
  n_grid = 0
endif

end function wavefront_shape

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_k0 (wf) result (k0)
!
! Routine to return the central wavenumber twopi / wavelength [rad/m].
!
! Input:
!   wf    -- wavefront_struct: Wavefront.
!
! Output:
!   k0    -- real(rp): Central wavenumber [rad/m].
!-

function wavefront_k0 (wf) result (k0)

type (wavefront_struct) wf
real(rp) k0

!

k0 = twopi / wf%wavelength

end function wavefront_k0

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_photon_energy (wf) result (e_photon)
!
! Routine to return the central photon energy h_bar * c * k0 [eV].
!
! Input:
!   wf        -- wavefront_struct: Wavefront.
!
! Output:
!   e_photon  -- real(rp): Central photon energy [eV].
!-

function wavefront_photon_energy (wf) result (e_photon)

type (wavefront_struct) wf
real(rp) e_photon

!

e_photon = wavefront_k0(wf) * h_bar_planck * c_light

end function wavefront_photon_energy

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_coord_vec (n_pt, d_pt) result (vec)
!
! Routine to return the coordinate vector of a zero centered uniform grid:
! n_pt points spaced d_pt apart running from -(n_pt-1)*d_pt/2 to +(n_pt-1)*d_pt/2.
!
! wavefront_xvec, wavefront_yvec and wavefront_zvec are this applied to the three axes.
!
! Input:
!   n_pt      -- integer: Number of points.
!   d_pt      -- real(rp): Point spacing.
!
! Output:
!   vec(n_pt) -- real(rp), allocatable: Coordinates.
!-

function wavefront_coord_vec (n_pt, d_pt) result (vec)

integer n_pt, i
real(rp) d_pt
real(rp), allocatable :: vec(:)

!

allocate (vec(n_pt))
do i = 1, n_pt
  vec(i) = (i - 1 - (n_pt - 1) / 2.0_rp) * d_pt
enddo

end function wavefront_coord_vec

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_xvec (wf) result (vec)
!
! Routine to return the x coordinates of the grid points [m].
!
! Input:
!   wf        -- wavefront_struct: Wavefront.
!
! Output:
!   vec(nx)   -- real(rp), allocatable: x coordinates [m].
!-

function wavefront_xvec (wf) result (vec)

type (wavefront_struct) wf
real(rp), allocatable :: vec(:)
integer n_grid(3)

!

n_grid = wavefront_shape(wf)
vec = wavefront_coord_vec(n_grid(1), wf%dx)

end function wavefront_xvec

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_yvec (wf) result (vec)
!
! Routine to return the y coordinates of the grid points [m].
!
! Input:
!   wf        -- wavefront_struct: Wavefront.
!
! Output:
!   vec(ny)   -- real(rp), allocatable: y coordinates [m].
!-

function wavefront_yvec (wf) result (vec)

type (wavefront_struct) wf
real(rp), allocatable :: vec(:)
integer n_grid(3)

!

n_grid = wavefront_shape(wf)
vec = wavefront_coord_vec(n_grid(2), wf%dy)

end function wavefront_yvec

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_zvec (wf) result (vec)
!
! Routine to return the z coordinates of the slices [m].
!
! Input:
!   wf        -- wavefront_struct: Wavefront.
!
! Output:
!   vec(nz)   -- real(rp), allocatable: z coordinates [m].
!-

function wavefront_zvec (wf) result (vec)

type (wavefront_struct) wf
real(rp), allocatable :: vec(:)
integer n_grid(3)

!

n_grid = wavefront_shape(wf)
vec = wavefront_coord_vec(n_grid(3), wf%dz)

end function wavefront_zvec

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_fft_wavenumber_vec (n_pt, d_pt) result (k_vec)
!
! Routine to return the transverse wavenumbers in FFT storage order. That is, in the order
! the values come out of a forward transform: zero frequency first, ascending positive
! frequencies, then the negative frequencies ascending to just below zero. This is
! twopi * numpy.fft.fftfreq(n_pt, d = d_pt). No fftshift is applied, since a propagation
! kernel built from these is used between a forward and a backward transform.
!
! Input:
!   n_pt        -- integer: Number of points.
!   d_pt        -- real(rp): Point spacing [m].
!
! Output:
!   k_vec(n_pt) -- real(rp), allocatable: Wavenumbers [rad/m].
!-

function wavefront_fft_wavenumber_vec (n_pt, d_pt) result (k_vec)

integer n_pt, i, m
real(rp) d_pt
real(rp), allocatable :: k_vec(:)

!

allocate (k_vec(n_pt))

do i = 1, n_pt
  m = i - 1
  if (2 * m > n_pt - 1) m = m - n_pt   ! The second half of the array holds negative wavenumbers.
  k_vec(i) = twopi * m / (n_pt * d_pt)
enddo

end function wavefront_fft_wavenumber_vec

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_intensity (wf) result (intensity)
!
! Routine to return the field intensity c * eps_0 / 2 * (|Ex|^2 + |Ey|^2) [W/m^2].
!
! Input:
!   wf                   -- wavefront_struct: Wavefront.
!
! Output:
!   intensity(nx,ny,nz)  -- real(rp), allocatable: Intensity.
!-

function wavefront_intensity (wf) result (intensity)

type (wavefront_struct) wf
real(rp), allocatable :: intensity(:,:,:)
integer n_grid(3)

!

n_grid = wavefront_shape(wf)
allocate (intensity(n_grid(1), n_grid(2), n_grid(3)))
intensity = 0

if (allocated(wf%Ex)) intensity = intensity + abs(wf%Ex)**2
if (allocated(wf%Ey)) intensity = intensity + abs(wf%Ey)**2

intensity = c_light * eps_0_vac / 2 * intensity

end function wavefront_intensity

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_energy (wf) result (energy)
!
! Routine to return the total field energy eps_0/2 * Integral |E|^2 dx dy dz [J].
!
! Note: The energy is conserved exactly by wavefront_drift. The propagation kernel has unit
! modulus everywhere, and the transform pair is unitary up to the 1/(nx*ny) factor. So
! the energy is a useful invariant to watch.
!
! Input:
!   wf      -- wavefront_struct: Wavefront.
!
! Output:
!   energy  -- real(rp): Total field energy [J].
!-

function wavefront_energy (wf) result (energy)

type (wavefront_struct) wf
real(rp) energy

!

energy = sum(wavefront_intensity(wf)) / c_light * wf%dx * wf%dy * wf%dz

end function wavefront_energy

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_fluence (wf) result (fluence)
!
! Routine to return the transverse fluence eps_0/2 * Integral |E(x,y,z)|^2 dz [J/m^2].
!
! Input:
!   wf             -- wavefront_struct: Wavefront.
!
! Output:
!   fluence(nx,ny) -- real(rp), allocatable: Fluence.
!-

function wavefront_fluence (wf) result (fluence)

type (wavefront_struct) wf
real(rp), allocatable :: fluence(:,:)

!

fluence = sum(wavefront_intensity(wf), dim = 3) / c_light * wf%dz

end function wavefront_fluence

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_power (wf) result (power)
!
! Routine to return the longitudinal power profile Integral I(x,y,z) dx dy [W].
!
! Input:
!   wf         -- wavefront_struct: Wavefront.
!
! Output:
!   power(nz)  -- real(rp), allocatable: Power per slice.
!-

function wavefront_power (wf) result (power)

type (wavefront_struct) wf
real(rp), allocatable :: power(:)

!

power = sum(sum(wavefront_intensity(wf), dim = 1), dim = 1) * wf%dx * wf%dy

end function wavefront_power

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_transverse_moments (wf, mean_x, mean_y, sigma_x, sigma_y)
!
! Routine to return the intensity weighted transverse centroid and rms size of a
! wavefront, summed over all slices. Mirrors the Python class's mean_x, mean_y, sigma_x
! and sigma_y.
!
! These give an independent analytic check on wavefront_drift. A Gaussian of waist size
! sigma0, drifted a distance z from its waist, has
!   sigma(z) = sqrt(sigma0^2 + (z * wavelength / (fourpi * sigma0))^2)
!
! Input:
!   wf        -- wavefront_struct: Wavefront.
!
! Output:
!   mean_x    -- real(rp): Intensity weighted <x> [m].
!   mean_y    -- real(rp): Intensity weighted <y> [m].
!   sigma_x   -- real(rp): Intensity weighted rms x [m].
!   sigma_y   -- real(rp): Intensity weighted rms y [m].
!-

subroutine wavefront_transverse_moments (wf, mean_x, mean_y, sigma_x, sigma_y)

type (wavefront_struct) wf
real(rp) mean_x, mean_y, sigma_x, sigma_y
real(rp) wt_sum
real(rp), allocatable :: intensity(:,:,:), wt_x(:), wt_y(:), xvec(:), yvec(:)

!

intensity = wavefront_intensity(wf)
wt_x = sum(sum(intensity, dim = 3), dim = 2)
wt_y = sum(sum(intensity, dim = 3), dim = 1)
xvec = wavefront_xvec(wf)
yvec = wavefront_yvec(wf)

wt_sum = sum(wt_x)
if (wt_sum <= 0) then
  mean_x = 0; mean_y = 0; sigma_x = 0; sigma_y = 0
  return
endif

mean_x = sum(wt_x * xvec) / wt_sum
mean_y = sum(wt_y * yvec) / wt_sum

! Two pass variance rather than <x^2> - <x>^2. The one pass form differences two numbers
! of similar size when the centroid is large against the width. That is normal for an
! off-axis beam, and the one pass form then loses most of its significant figures.
! mean_variance_calc in openPMD-beamphysics uses the same two pass form.

sigma_x = sqrt(sum(wt_x * (xvec - mean_x)**2) / wt_sum)
sigma_y = sqrt(sum(wt_y * (yvec - mean_y)**2) / wt_sum)

end subroutine wavefront_transverse_moments

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_drift (wf, z_drift, err_flag, curvature)
!
! Routine to propagate a wavefront a distance z_drift through free space.
!
! Each slice is propagated independently with the exact transfer function of the paraxial
! wave equation, applied in transverse k-space:
!
!   E(kx, ky) -> E(kx, ky) * exp(-i * z_drift * wavelength * (kx^2 + ky^2) / fourpi)
!
! The exponent is -i * z_drift * (kx^2 + ky^2) / (2 * k0) written in terms of the
! wavelength, since wavelength / fourpi = 1 / (2 * k0). This is spectrally exact and
! unconditionally stable. It is not a finite difference step, and subdividing z_drift
! gains no accuracy. It is the same kernel as Genesis's FieldSolverFFT and as
! drift_wavefront_basic in openPMD-beamphysics.
!
! The transverse boundary is periodic, since a discrete Fourier transform imposes that and
! no apodisation is applied. Field reaching the transverse edge of the grid wraps around
! and re-enters. Over a long drift the grid must be large enough to contain the field, or
! be padded until it is.
!
! s_position advances by z_drift. The mesh's own z axis is the intra-pulse coordinate and
! is co-moving. So the kernel leaves it alone, and s_position is the only place a
! propagation along the beamline is recorded. openPMD-beamphysics splits the two the same
! way, and its drift advances s_position too. Operations that do not move the pulse down
! the line (crop, pad, a lens phase) leave s_position alone.
!
! Note: z_drift is the optical path. Through a geometry break the light goes the chord
! while the lattice position advances by the arc. So the FEL walk, which knows the
! lattice, sets s_position itself after each piece (fel_track_line_mod).
!
! Input:
!   wf          -- wavefront_struct: Wavefront to propagate.
!   z_drift     -- real(rp): Drift distance [m]. May be negative.
!   curvature   -- real(rp), optional: Accepted only as zero. Default is 0. A nonzero
!                    curvature selects quadratic phase rescaling with an expanding grid,
!                    which is not implemented here. Passing one is an error.
!
! Output:
!   wf          -- wavefront_struct: Propagated wavefront, s_position advanced.
!   err_flag    -- logical, optional: Set True on error, False otherwise.
!-

subroutine wavefront_drift (wf, z_drift, err_flag, curvature)

type (wavefront_struct), target :: wf
real(rp) z_drift
real(rp), optional :: curvature
logical, optional :: err_flag
logical err

integer n_grid(3), nx, ny, nz, iz, i_pol
complex(wf_rp), allocatable :: kernel(:,:)
complex(wf_rp), pointer :: field(:,:,:)
logical any_err
character(*), parameter :: r_name = 'wavefront_drift'

!

if (present(err_flag)) err_flag = .true.

call wavefront_check (wf, err);  if (err) return

if (present(curvature)) then
  if (curvature /= 0) then
    call out_io (s_error$, r_name, 'NONZERO CURVATURE IS NOT IMPLEMENTED BY THIS PROPAGATOR.', &
                 'Curvature corrected propagation with grid rescaling, as in', &
                 'drift_wavefront_advanced of openPMD-beamphysics, is a separate propagator.')
    return
  endif
endif

! Note: err_flag stays true until the transforms have finished. Clearing it early and then
! returning from the loop below on a failed transform would report success on failure.

if (z_drift == 0) then
  if (present(err_flag)) err_flag = .false.
  return
endif

n_grid = wavefront_shape(wf)
nx = n_grid(1); ny = n_grid(2); nz = n_grid(3)

kernel = wavefront_drift_kernel(wf, z_drift)

! Apply slice by slice. Both polarisation components see the same kernel, since the
! paraxial operator does not couple Ex and Ey in free space.
!
! Parallel over slices, as in the FEL field solve (fel_field_step). Every slice touches
! only its own plane, the kernel is read-only, and the FFTW plan cache is threadprivate.
! So the arithmetic per slice does not depend on which thread runs it or on how many
! there are. A time-dependent line makes one of these calls per field-free element.
! At production slice counts a serial loop here is the largest serial block of the walk
! (doc/performance.md perf-thread-scaling).
!
! wavefront_fft2_plan_threads warms every thread's plans from the serial context first.
! FFTW's planner is not reentrant, and creating a plan inside this loop would race with
! the transforms.

call wavefront_fft2_plan_threads (nx, ny, err);  if (err) return

do i_pol = 1, 2
  if (i_pol == 1) then
    if (.not. allocated(wf%Ex)) cycle
    field => wf%Ex
  else
    if (.not. allocated(wf%Ey)) cycle
    field => wf%Ey
  endif

  any_err = .false.
  !$OMP parallel do private(err) reduction(.or.: any_err)
  do iz = 1, nz
    call wavefront_fft2 (field(:,:,iz), wf_fft_forward$, err)
    if (.not. err) then
      field(:,:,iz) = kernel * field(:,:,iz)
      call wavefront_fft2 (field(:,:,iz), wf_fft_backward$, err)
      if (.not. err) field(:,:,iz) = field(:,:,iz) / real(nx * ny, wf_rp)
    endif
    any_err = any_err .or. err
  enddo
  !$OMP end parallel do
  if (any_err) return
enddo

wf%s_position = wf%s_position + z_drift

if (present(err_flag)) err_flag = .false.

end subroutine wavefront_drift

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_drift_kernel (wf, z_drift) result (kernel)
!
! Routine to return the transverse k-space propagation kernel for a drift of z_drift:
!
!   exp(-i * z_drift * wavelength * (kx^2 + ky^2) / fourpi)
!
! in FFT storage order, so it can be applied between a forward and a backward transform
! with no fftshift.
!
! The arithmetic is written in the same association order as drift_wavefront_basic in
! openPMD-beamphysics. So the two implementations differ by round-off and not by grouping.
!
! Input:
!   wf              -- wavefront_struct: Supplies nx, ny, dx, dy and wavelength.
!   z_drift         -- real(rp): Drift distance [m].
!
! Output:
!   kernel(nx,ny)   -- complex(wf_rp), allocatable: Propagation kernel.
!-

function wavefront_drift_kernel (wf, z_drift) result (kernel)

type (wavefront_struct) wf
real(rp) z_drift
complex(wf_rp), allocatable :: kernel(:,:)
integer n_grid(3), nx, ny, ix, iy
real(rp), allocatable :: kx_vec(:), ky_vec(:)
real(rp) phase, k2_scale

!

n_grid = wavefront_shape(wf)
nx = n_grid(1); ny = n_grid(2)

kx_vec = wavefront_fft_wavenumber_vec(nx, wf%dx)
ky_vec = wavefront_fft_wavenumber_vec(ny, wf%dy)
k2_scale = -z_drift * wf%wavelength

allocate (kernel(nx, ny))
do iy = 1, ny
  do ix = 1, nx
    phase = k2_scale * (kx_vec(ix)**2 + ky_vec(iy)**2) / (4 * pi)
    kernel(ix, iy) = cmplx(cos(phase), sin(phase), wf_rp)
  enddo
enddo

end function wavefront_drift_kernel

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_drift_reference (wf, z_drift, err_flag)
!
! Routine to propagate a wavefront as wavefront_drift does, sharing no code with it.
!
! This is a validation reference and not a propagator for use. It costs
! O(nx*ny*(nx+ny)) per slice against the FFT's O(nx*ny*log(nx*ny)). It builds its own
! propagation kernel and performs its own transform. So it is independent of
! wavefront_drift_kernel and of every FFTW convention: the dimension order handed to
! fftw_plan_dft_2d, the sign of each direction, the normalization factor, and the
! alignment rules for executing a cached plan on a differently allocated array.
!
! The comparison against openPMD-beamphysics goes through the Genesis4 file format, which
! requires nx = ny and dx = dy. On a square grid with equal spacings the propagation
! kernel is symmetric under interchanging the two transverse axes. So that comparison
! cannot see a transposed transform, nor an interchanged dx and dy, however asymmetric
! the test field is. Both errors pass the Python comparison at round-off. Only a
! rectangular grid with unequal spacings can see them, and only this routine supplies a
! reference on one.
!
! The kernel below is written in a different but equivalent form from
! wavefront_drift_kernel. The wavenumbers are formed inline rather than through
! wavefront_fft_wavenumber_vec, and the coefficient is 1/(2*k0) rather than
! wavelength/fourpi. So the comparison checks the arithmetic as well as the transform.
! The two differ by round-off in the phase as well as in the transform.
!
! Note: A direct sum accumulates round-off as O(n) rather than the FFT's O(log n), so
! this routine is the less accurate of the two. Expect agreement near 1e-14 rather than
! 1e-16. The residual is not an error in wavefront_drift.
!
! Input:
!   wf          -- wavefront_struct: Wavefront to propagate.
!   z_drift     -- real(rp): Drift distance [m].
!
! Output:
!   wf          -- wavefront_struct: Propagated wavefront.
!   err_flag    -- logical, optional: Set True on error, False otherwise.
!-

subroutine wavefront_drift_reference (wf, z_drift, err_flag)

type (wavefront_struct), target :: wf
real(rp) z_drift
logical, optional :: err_flag
logical err

integer n_grid(3), nx, ny, nz, ix, iy, iz, i_pol, mx, my
complex(wf_rp), allocatable :: kernel(:,:), slice(:,:)
complex(wf_rp), pointer :: field(:,:,:)
real(rp) k0, kx, ky, phase

!

if (present(err_flag)) err_flag = .true.
call wavefront_check (wf, err);  if (err) return

! As in wavefront_drift, err_flag is cleared only once the work is done. Nothing below can
! fail, since the direct sum has no fallible call. Clearing it early would report success
! on failure once one is added.

if (z_drift == 0) then
  if (present(err_flag)) err_flag = .false.
  return
endif

n_grid = wavefront_shape(wf)
nx = n_grid(1); ny = n_grid(2); nz = n_grid(3)

! The kernel is built here rather than shared (see the header).

k0 = twopi / wf%wavelength
allocate (kernel(nx, ny), slice(nx, ny))

do iy = 1, ny
  my = iy - 1
  if (my > (ny - 1) / 2) my = my - ny
  ky = my * twopi / (ny * wf%dy)

  do ix = 1, nx
    mx = ix - 1
    if (mx > (nx - 1) / 2) mx = mx - nx
    kx = mx * twopi / (nx * wf%dx)

    phase = -z_drift * (kx * kx + ky * ky) / (2 * k0)
    kernel(ix, iy) = cmplx(cos(phase), sin(phase), wf_rp)
  enddo
enddo

do i_pol = 1, 2
  if (i_pol == 1) then
    if (.not. allocated(wf%Ex)) cycle
    field => wf%Ex
  else
    if (.not. allocated(wf%Ey)) cycle
    field => wf%Ey
  endif

  do iz = 1, nz
    ! Forward: transform along x, then along y.
    do iy = 1, ny
      slice(:,iy) = wavefront_dft_1d(field(:,iy,iz), wf_fft_forward$)
    enddo
    do ix = 1, nx
      slice(ix,:) = wavefront_dft_1d(slice(ix,:), wf_fft_forward$)
    enddo

    slice = kernel * slice

    ! Backward, same two passes with the opposite sign, then the 1/(nx*ny) normalization.
    do iy = 1, ny
      slice(:,iy) = wavefront_dft_1d(slice(:,iy), wf_fft_backward$)
    enddo
    do ix = 1, nx
      slice(ix,:) = wavefront_dft_1d(slice(ix,:), wf_fft_backward$)
    enddo

    field(:,:,iz) = slice / real(nx * ny, wf_rp)
  enddo
enddo

if (present(err_flag)) err_flag = .false.

end subroutine wavefront_drift_reference

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_dft_1d (dat, direction) result (out)
!
! Routine to return the unnormalized discrete Fourier transform of a vector, computed as a
! direct sum:
!
!   out(k) = Sum_j dat(j) * exp(direction * i * twopi * (j-1) * (k-1) / n)
!
! direction = wf_fft_forward$ gives the negative exponent, matching FFTW and numpy.fft.
! Used only by wavefront_drift_reference. See the discussion there.
!
! Input:
!   dat(:)      -- complex(wf_rp): Data to transform.
!   direction   -- integer: wf_fft_forward$ or wf_fft_backward$.
!
! Output:
!   out(size(dat)) -- complex(wf_rp), allocatable: Transform.
!-

function wavefront_dft_1d (dat, direction) result (out)

complex(wf_rp) dat(:)
complex(wf_rp), allocatable :: out(:)
integer direction, n_pt, j, k
integer(8) n8
real(rp) angle

!

n_pt = size(dat)
n8 = n_pt
allocate (out(n_pt))

! The index product is reduced modulo n before it becomes an angle. That keeps the angle
! inside one turn, so cos and sin need no argument reduction. The 8-byte product also
! cannot overflow at large n as a default integer would.

do k = 1, n_pt
  out(k) = 0
  do j = 1, n_pt
    angle = direction * twopi * modulo(int(j-1, 8) * int(k-1, 8), n8) / n_pt
    out(k) = out(k) + dat(j) * cmplx(cos(angle), sin(angle), wf_rp)
  enddo
enddo

end function wavefront_dft_1d

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_fft2 (dat, direction, err_flag)
!
! Routine to apply an unnormalized in-place 2D complex FFT to dat using a cached FFTW plan.
!
! Neither direction applies a normalization factor. So a forward followed by a backward
! transform multiplies the data by nx*ny. This matches FFTW and numpy.fft.fft2. numpy's
! ifft2 differs from wf_fft_backward$ only by that factor.
!
! Plans are created once per grid size and cached per thread, following
! bmad/space_charge/fft_interface_mod.f90. The transform is executed on an internal
! fftw_alloc_complex buffer rather than on dat. A plan may only be executed on an array
! with the alignment it was created with. dat is copied in and out.
!
! Note: The plan cache and its work buffer are threadprivate, so concurrent calls from a
! parallel loop over slices each use their own. A caller running transforms in a parallel
! loop must first warm every thread's cache with wavefront_fft2_plan_threads. A lazy
! first-call plan would run the planner concurrently with other threads' transforms,
! which FFTW does not promise to be safe. Also see wavefront_fft2_plan.
!
! Input:
!   dat(:,:)    -- complex(wf_rp): Data to transform.
!   direction   -- integer: wf_fft_forward$ or wf_fft_backward$.
!
! Output:
!   dat(:,:)    -- complex(wf_rp): Transformed data.
!   err_flag    -- logical, optional: Set True on error, False otherwise.
!-

subroutine wavefront_fft2 (dat, direction, err_flag)

include 'fftw3.f03'

complex(wf_rp) dat(:,:)
integer direction
logical, optional :: err_flag
logical cache_ok
integer nx, ny
character(*), parameter :: r_name = 'wavefront_fft2'

!

if (present(err_flag)) err_flag = .true.

if (direction /= wf_fft_forward$ .and. direction /= wf_fft_backward$) then
  call out_io (s_error$, r_name, 'BAD TRANSFORM DIRECTION: \i0\ ', i_array = [direction])
  return
endif

nx = size(dat, 1)
ny = size(dat, 2)

call wavefront_fft2_plan (nx, ny, cache_ok)

if (.not. cache_ok) then
  call out_io (s_error$, r_name, &
        'FFTW COULD NOT ALLOCATE A WORK BUFFER OR CREATE A PLAN FOR A \i0\ BY \i0\ TRANSFORM.', &
        i_array = [nx, ny])
  return
endif

wf_buf = dat

! Passing the same array as both in and out is FFTW's documented way to ask for an
! in-place transform, as bmad/space_charge/fft_interface_mod.f90 does. gfortran -Wall
! warns about it, because the fftw3.f03 interface declares both dummies intent(out).
! Aliasing those would be nonconforming for a normal Fortran procedure. This one is
! bind(C) with assumed-size dummies, so the arguments are passed by address with no copy
! in or out, and the aliasing is safe. Bmad's own build does not enable that warning.

if (direction == wf_fft_forward$) then
  call fftw_execute_dft (wf_plan_fwd, wf_buf, wf_buf)
else
  call fftw_execute_dft (wf_plan_bwd, wf_buf, wf_buf)
endif

dat = wf_buf

if (present(err_flag)) err_flag = .false.

end subroutine wavefront_fft2

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_fft2_plan (nx, ny, cache_ok)
!
! Routine to fill the calling thread's plan cache for an nx by ny transform. Nothing is
! done when the cache already holds this size. The critical section guards FFTW's rule
! that only fftw_execute is reentrant. So the planner is globally serialised, even though
! every thread builds into its own threadprivate cache.
!
! The same rule is why wavefront_fft2_plan_threads exists. A thread that plans lazily on
! its first transform does so while other threads are executing transforms. Planner
! activity concurrent with execution is outside FFTW's thread-safety promise. It makes
! the planes of the last threads to plan differ at the ulp level from run to run. A
! caller about to run transforms in a parallel loop must warm every thread's cache first.
!
! A failure is recorded in wf_cache_nx and reported by the caller, rather than signalled
! from inside the section. Branching out of a critical construct is not allowed.
!
! Input:
!   nx, ny      -- integer: Transform size.
!
! Output:
!   cache_ok    -- logical: True if this thread's cache now holds plans for (nx, ny).
!-

subroutine wavefront_fft2_plan (nx, ny, cache_ok)

include 'fftw3.f03'

integer nx, ny
logical cache_ok

!

!$OMP CRITICAL (wavefront_fft_plan_lock)

if (nx /= wf_cache_nx .or. ny /= wf_cache_ny) then
  if (C_ASSOCIATED(wf_plan_fwd)) call fftw_destroy_plan(wf_plan_fwd)
  if (C_ASSOCIATED(wf_plan_bwd)) call fftw_destroy_plan(wf_plan_bwd)
  if (C_ASSOCIATED(wf_buf_cptr)) call fftw_free(wf_buf_cptr)
  wf_plan_fwd = C_NULL_PTR
  wf_plan_bwd = C_NULL_PTR
  wf_buf_cptr = C_NULL_PTR
  wf_buf => null()

  ! Mark the cache invalid first, so a failure below leaves it invalid rather than half
  ! built. It is set to the real size only once every piece is obtained.
  wf_cache_nx = 0; wf_cache_ny = 0

  ! fftw_alloc_complex and fftw_plan_dft_2d both return null on failure. A null plan
  ! reaching fftw_execute_dft crashes rather than giving a diagnosable error.
  wf_buf_cptr = fftw_alloc_complex(int(nx, C_SIZE_T) * int(ny, C_SIZE_T))

  if (C_ASSOCIATED(wf_buf_cptr)) then
    call C_F_POINTER (wf_buf_cptr, wf_buf, [nx, ny])

    ! Note the reversed dimension order. fftw_plan_dft_2d takes the slowest varying
    ! dimension first, and in a Fortran (nx, ny) array that is ny. Getting this backwards
    ! is invisible on a square grid, which is why wavefront_drift_reference exists.
    !
    ! FFTW_ESTIMATE is used and not FFTW_MEASURE. FFTW_MEASURE picks the algorithm by
    ! timing candidate transforms, so two runs of the same binary can pick differently. Every
    ! transform after that differs at the ulp level, enough to change a run's field
    ! diagnostics byte for byte. ESTIMATE's choice is a pure function of the problem. So
    ! results reproduce from run to run, from one thread count to another, and against the
    ! banked reference files. Reproducibility matters more here than FFT speed.
    wf_plan_fwd = fftw_plan_dft_2d(ny, nx, wf_buf, wf_buf, FFTW_FORWARD,  FFTW_ESTIMATE)
    wf_plan_bwd = fftw_plan_dft_2d(ny, nx, wf_buf, wf_buf, FFTW_BACKWARD, FFTW_ESTIMATE)

    if (C_ASSOCIATED(wf_plan_fwd) .and. C_ASSOCIATED(wf_plan_bwd)) then
      wf_cache_nx = nx; wf_cache_ny = ny
    endif
  endif
endif

cache_ok = (wf_cache_nx == nx .and. wf_cache_ny == ny)

!$OMP END CRITICAL (wavefront_fft_plan_lock)

end subroutine wavefront_fft2_plan

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_fft2_plan_threads (nx, ny, err_flag)
!
! Routine to fill every OpenMP thread's plan cache for an nx by ny transform. No planner
! then runs concurrently with transform execution, since only fftw_execute is reentrant
! (also see wavefront_fft2_plan). Call it serially, before any parallel loop that executes
! transforms of this size. Nothing is done when every thread already holds this size.
!
! Input:
!   nx, ny      -- integer: Transform size.
!
! Output:
!   err_flag    -- logical: Set True if any thread's planning failed, False otherwise.
!-

subroutine wavefront_fft2_plan_threads (nx, ny, err_flag)

integer nx, ny
logical err_flag
logical any_bad, cache_ok
character(*), parameter :: r_name = 'wavefront_fft2_plan_threads'

!

any_bad = .false.

!$OMP PARALLEL private(cache_ok) reduction(.or.: any_bad)
call wavefront_fft2_plan (nx, ny, cache_ok)
any_bad = any_bad .or. .not. cache_ok
!$OMP END PARALLEL

if (any_bad) then
  call out_io (s_error$, r_name, &
        'FFTW COULD NOT ALLOCATE A WORK BUFFER OR CREATE A PLAN FOR A \i0\ BY \i0\ TRANSFORM.', &
        i_array = [nx, ny])
endif

err_flag = any_bad

end subroutine wavefront_fft2_plan_threads

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_fft_free ()
!
! Routine to destroy the cached FFTW plans and free the work buffer of the calling
! thread. The cache is threadprivate. So a serial call after parallel work leaves the
! worker threads' caches allocated until program end. Freeing those would need a call
! from inside a parallel region. Not needed for correctness. It makes a single-threaded
! leak check clean.
!-

subroutine wavefront_fft_free ()

include 'fftw3.f03'

!

!$OMP CRITICAL (wavefront_fft_plan_lock)
if (C_ASSOCIATED(wf_plan_fwd)) call fftw_destroy_plan(wf_plan_fwd)
if (C_ASSOCIATED(wf_plan_bwd)) call fftw_destroy_plan(wf_plan_bwd)
if (C_ASSOCIATED(wf_buf_cptr)) call fftw_free(wf_buf_cptr)
wf_plan_fwd = C_NULL_PTR
wf_plan_bwd = C_NULL_PTR
wf_buf_cptr = C_NULL_PTR
wf_buf => null()
wf_cache_nx = 0; wf_cache_ny = 0
!$OMP END CRITICAL (wavefront_fft_plan_lock)

end subroutine wavefront_fft_free

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_params_of_plane (plane, dx, wavelength, dz_slice, pms, with_angles, err_flag, &
!                                                                        theta_cut, power_inside)
!
! Routine to compute the wavefront_params of one transverse field plane (one slice).
!
! Spatial moments come from the intensity |E|^2 in one grid pass. The theta moments use
! the spectral representation, with theta = k_perp/ks. <theta^2> comes from |FFT(E)|^2.
! The cross moment is <x theta>_sym = Im INT conj(E) x dE/dx / (ks INT |E|^2), with the
! derivative evaluated spectrally. That takes one forward and two inverse FFTs. So angle
! moments are computed only where needed (with_angles, recorded by angle_moments_valid).
! The grid k ordering matches the field solver's kernel (fftshift index mapping).
!
! Input:
!   plane(:,:)    -- complex(wf_rp): One field slice [V/m].
!   dx            -- real(rp): Transverse grid spacing [m] (square grid).
!   wavelength    -- real(rp): Radiation wavelength [m].
!   dz_slice      -- real(rp): Slice spacing [m] (energy = power * dz_slice/c).
!   with_angles   -- logical: Fill the theta rows of sigma? Costs 3 FFTs.
!   theta_cut     -- real(rp), optional: Angular radius [rad] for power_inside.
!
! Output:
!   pms           -- wavefront_params_struct: The moments. pms%s is *not* set here.
!   err_flag      -- logical: Set True on error (FFT failure), False otherwise.
!   power_inside  -- real(rp), optional: Power [W] within theta_cut of the axis. NaN when
!                      with_angles is False or the plane carries no power.
!-

subroutine wavefront_params_of_plane (plane, dx, wavelength, dz_slice, pms, with_angles, err_flag, &
                                                                        theta_cut, power_inside)

complex(wf_rp) plane(:,:)
real(rp) dx, wavelength, dz_slice
type (wavefront_params_struct) pms
real(rp) nan
logical with_angles, err_flag
real(rp), optional :: theta_cut, power_inside
real(rp) s_in, theta2, cut2

complex(wf_rp), allocatable :: ft(:,:), dxe(:,:), dye(:,:)
real(rp) wt, wsum, x, y, ks, dk, shift, kx, ky, scl
real(rp) sx, sy, sxx, syy, stx, sty, stxx, styy, sxtx, syty
integer nx, ny, ix, iy, ic
logical err

!

err_flag = .true.
pms = wavefront_params_struct()
nx = size(plane, 1);  ny = size(plane, 2)
ks = twopi / wavelength

! Intensity moments in one separable pass. The inner loop only accumulates row and
! column intensity sums, which vectorizes. The 1-D moments follow. Coordinates are
! grid-centered, matching the solver. This runs per slice per record, so its cost is the
! overhead of the stats file. Keep it lean.

shift = -0.5_rp * (nx - 1)
block
  real(rp) colsum(nx), rowsum(ny), rs
  colsum = 0
  do iy = 1, ny
    rs = 0
    do ix = 1, nx
      wt = real(plane(ix,iy), rp)**2 + aimag(plane(ix,iy))**2
      rs = rs + wt
      colsum(ix) = colsum(ix) + wt
    enddo
    rowsum(iy) = rs
  enddo
  wsum = sum(rowsum)
  sx = 0;  sxx = 0
  do ix = 1, nx
    x = (ix - 1 + shift) * dx
    sx = sx + colsum(ix) * x;  sxx = sxx + colsum(ix) * x**2
  enddo
  sy = 0;  syy = 0
  do iy = 1, ny
    y = (iy - 1 + shift) * dx
    sy = sy + rowsum(iy) * y;  syy = syy + rowsum(iy) * y**2
  enddo
end block

scl = dx**2 / (2 * (mu_0_vac * c_light))       ! |E|^2 sum -> power [W], the diag convention.
pms%power = wsum * scl
pms%energy = pms%power * dz_slice / c_light
ic = nx/2 + 1
pms%on_axis_intensity = (real(plane(ic,ic), rp)**2 + aimag(plane(ic,ic))**2) / (2 * (mu_0_vac * c_light))

if (wsum > 0) then
  pms%centroid(1) = sx / wsum
  pms%centroid(3) = sy / wsum
  pms%sigma(1,1) = sxx / wsum - pms%centroid(1)**2
  pms%sigma(3,3) = syy / wsum - pms%centroid(3)**2
endif

! Angle moments not computed are set to NaN rather than zero. Zero is a legal value for a
! centroid and for a covariance entry, so a consumer could not tell it from an answer.
! angle_moments_valid says which case this is, and the NaN says it in the numbers too.
! A NaN also propagates through whatever a consumer does with it. Bmad's real_garbage$
! sentinel would plot as a number.

if (.not. with_angles .or. wsum <= 0) then
  nan = ieee_value(1.0_rp, ieee_quiet_nan)
  pms%centroid(2) = nan
  pms%centroid(4) = nan
  pms%sigma(2,:) = nan;  pms%sigma(:,2) = nan
  pms%sigma(4,:) = nan;  pms%sigma(:,4) = nan
  pms%emit_x = nan;      pms%emit_y = nan
  if (present(power_inside)) power_inside = nan
  err_flag = .false.
  return
endif

! Spectral side: FFT once for |ft|^2 (theta first and second moments), inverse-transform
! i*k*ft for the spatial derivatives feeding the cross moments.

allocate (ft(nx,ny), dxe(nx,ny), dye(nx,ny))
ft = plane
call wavefront_fft2 (ft, wf_fft_forward$, err);  if (err) return

dk = twopi / (nx * dx)
stx = 0;  sty = 0;  stxx = 0;  styy = 0

! The power inside an angular radius is summed on this transform rather than on one of its
! own. It is the quantity doc/startup-noise.md splits on. Inside it lies the radiation
! that can couple to the mode. Outside it lies the wide-angle emission of the point
! beamlets, which is most of an unseeded run's power at the grids and loads in common use.

s_in = 0
cut2 = 0
if (present(theta_cut)) cut2 = (theta_cut * ks)**2
! The kernel's fftshift storage reduces to the standard FFT frequency order:
! stored index s (0-based) holds frequency s for s <= (n-1)/2, s-n above.

do iy = 1, ny
  ky = (mod((iy-1) + (ny-1)/2, ny) - (ny-1)/2) * dk
  do ix = 1, nx
    kx = (mod((ix-1) + (nx-1)/2, nx) - (nx-1)/2) * dk
    wt = real(ft(ix,iy), rp)**2 + aimag(ft(ix,iy))**2
    stx = stx + wt * kx;   sty = sty + wt * ky
    stxx = stxx + wt * kx**2;   styy = styy + wt * ky**2
    dxe(ix,iy) = cmplx(0.0_rp, kx, wf_rp) * ft(ix,iy)
    dye(ix,iy) = cmplx(0.0_rp, ky, wf_rp) * ft(ix,iy)
    theta2 = kx**2 + ky**2
    if (theta2 <= cut2) s_in = s_in + wt
  enddo
enddo
! Parseval: sum|ft|^2 = nx*ny * sum|E|^2 for this unnormalized transform.
scl = wsum * real(nx, rp) * real(ny, rp)
pms%centroid(2) = stx / (scl * ks)
pms%centroid(4) = sty / (scl * ks)
pms%sigma(2,2) = stxx / (scl * ks**2) - pms%centroid(2)**2
pms%sigma(4,4) = styy / (scl * ks**2) - pms%centroid(4)**2

! Parseval again: sum|ft|^2 over a region, divided by nx*ny, is the |E|^2 sum of that
! region. The same dx^2/(2 Z0) makes it a power.

if (present(power_inside)) power_inside = s_in / (real(nx, rp) * real(ny, rp)) &
                                          * dx**2 / (2 * (mu_0_vac * c_light))

call wavefront_fft2 (dxe, wf_fft_backward$, err);  if (err) return
call wavefront_fft2 (dye, wf_fft_backward$, err);  if (err) return
dxe = dxe / real(nx*ny, rp)                        ! dE/dx on the grid.
dye = dye / real(nx*ny, rp)

! <x theta_x>_sym = Im INT conj(E) x dE/dx / (ks INT|E|^2), then recentered.

sxtx = 0;  syty = 0
do iy = 1, ny
  y = (iy - 1 + shift) * dx
  do ix = 1, nx
    x = (ix - 1 + shift) * dx
    sxtx = sxtx + x * aimag(conjg(plane(ix,iy)) * dxe(ix,iy))
    syty = syty + y * aimag(conjg(plane(ix,iy)) * dye(ix,iy))
  enddo
enddo
pms%sigma(1,2) = sxtx / (wsum * ks) - pms%centroid(1) * pms%centroid(2)
pms%sigma(2,1) = pms%sigma(1,2)
pms%sigma(3,4) = syty / (wsum * ks) - pms%centroid(3) * pms%centroid(4)
pms%sigma(4,3) = pms%sigma(3,4)

pms%emit_x = sqrt(max(0.0_rp, pms%sigma(1,1) * pms%sigma(2,2) - pms%sigma(1,2)**2))
pms%emit_y = sqrt(max(0.0_rp, pms%sigma(3,3) * pms%sigma(4,4) - pms%sigma(3,4)**2))
pms%angle_moments_valid = .true.

err_flag = .false.

end subroutine wavefront_params_of_plane

end module wavefront_mod
