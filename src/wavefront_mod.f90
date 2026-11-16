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
use wavefront_wisdom_mod

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

! How the plans are made (wavefront_fft2_set_plan). wf_plan_wisdom$ makes the plan the
! wisdom text of wavefront_wisdom_mod records, and where that text gives none, the
! estimated one. wf_plan_estimate$ is FFTW_ESTIMATE, a pure function of the problem.
! wf_plan_measure$ is FFTW_MEASURE, which times candidates and so can choose differently
! in two runs. The first two give the same plan in every run and on every thread. That
! is what keeps a run's digits the same from run to run and from one thread count to
! another, since two plans for one transform differ in the last bit of every output.
!
! wf_plan_mode is shared and is set before any thread plans. wf_cache_mode and
! wf_cache_made belong to the calling thread's cache: the mode it was filled under, and
! the way its plans were in fact made, which is wf_plan_estimate$ under wf_plan_wisdom$
! where the wisdom gave no plan.

integer, parameter :: wf_plan_wisdom$ = 1, wf_plan_estimate$ = 2, wf_plan_measure$ = 3
integer, parameter :: wf_plan_wisdom_file$ = 4
character(12), parameter :: wf_plan_name(4) = [character(12) :: 'wisdom', 'estimate', 'measure', 'wisdom_file']

! What wavefront_fft2_make_wisdom found for a grid.

integer, parameter :: wf_wisdom_in_source$ = 1, wf_wisdom_in_file$ = 2, wf_wisdom_measured$ = 3

! The wisdom of a file, read once by wavefront_fft2_set_wisdom_file and held as a C
! string. Under wf_plan_wisdom$ a transform the text of wavefront_wisdom_mod gives no
! plan for is planned from this text, where it holds one, before the estimate. The text
! of the source is tried first, so a transform it holds is planned the same whatever
! file a machine carries. Shared, and set before any thread plans.

character(:, kind = C_CHAR), allocatable, private, save :: wf_file_wisdom

integer, private, save :: wf_plan_mode = wf_plan_wisdom$
integer, private, save :: wf_cache_mode = 0, wf_cache_made = 0
!$OMP threadprivate(wf_cache_mode, wf_cache_made)

type(C_PTR), private, save :: wf_plan_fwd = C_NULL_PTR
type(C_PTR), private, save :: wf_plan_bwd = C_NULL_PTR
type(C_PTR), private, save :: wf_buf_cptr = C_NULL_PTR
complex(wf_rp), private, pointer, contiguous, save :: wf_buf(:,:) => null()

! A second buffer of the same size and alignment, for a caller that holds two planes at
! once (wavefront_fft2_buffers). It is allocated at the first request and freed with
! the cache.

type(C_PTR), private, save :: wf_buf2_cptr = C_NULL_PTR
complex(wf_rp), private, pointer, contiguous, save :: wf_buf2(:,:) => null()
!$OMP threadprivate(wf_buf2_cptr, wf_buf2)
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
! (docs/performance.md perf-thread-scaling).
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
logical err
complex(wf_rp), pointer, contiguous :: buf(:,:)
character(*), parameter :: r_name = 'wavefront_fft2'

!

if (present(err_flag)) err_flag = .true.

if (direction /= wf_fft_forward$ .and. direction /= wf_fft_backward$) then
  call out_io (s_error$, r_name, 'BAD TRANSFORM DIRECTION: \i0\ ', i_array = [direction])
  return
endif

call wavefront_fft2_buffers (size(dat, 1), size(dat, 2), buf, err_flag = err)
if (err) return

! The transform runs on the thread's own buffer and not on the caller's array, so the
! plan always sees the alignment it was created with. A caller that would save the two
! copies works in the buffer itself (wavefront_fft2_buffers).

buf = dat
call wavefront_fft2_aligned (1, direction, err)
if (err) return
dat = buf

if (present(err_flag)) err_flag = .false.

end subroutine wavefront_fft2

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_fft2_buffers (nx, ny, buf, buf2, err_flag)
!
! Routine to return the calling thread's work buffers for an nx by ny transform, making
! the thread's plans if its cache does not hold this size. The buffers are allocated by
! fftw_alloc_complex, so wavefront_fft2_aligned may transform either of them where it
! lies, with no copy. A caller that transforms a plane, works on the result and
! transforms back does all of it in a buffer, and copies the plane in once and out once.
!
! The buffers are the thread's until its next call of wavefront_fft2, which uses the
! first, or of any routine that transforms through it. They hold whatever was left in
! them. A change of size frees both.
!
! Input:
!   nx, ny    -- integer: Transform size.
!
! Output:
!   buf       -- complex(wf_rp), pointer, contiguous (:,:): The first buffer, nx by ny.
!   buf2      -- complex(wf_rp), pointer, contiguous (:,:), optional: A second buffer of
!                  the same size.
!   err_flag  -- logical, optional: Set True if a buffer or a plan could not be made.
!-

subroutine wavefront_fft2_buffers (nx, ny, buf, buf2, err_flag)

include 'fftw3.f03'

integer nx, ny
complex(wf_rp), pointer, contiguous :: buf(:,:)
complex(wf_rp), pointer, contiguous, optional :: buf2(:,:)
logical, optional :: err_flag
logical cache_ok
character(*), parameter :: r_name = 'wavefront_fft2_buffers'

!

if (present(err_flag)) err_flag = .true.
buf => null()
if (present(buf2)) buf2 => null()

call wavefront_fft2_plan (nx, ny, cache_ok)

if (cache_ok .and. present(buf2) .and. .not. C_ASSOCIATED(wf_buf2_cptr)) then
  wf_buf2_cptr = fftw_alloc_complex(int(nx, C_SIZE_T) * int(ny, C_SIZE_T))
  if (C_ASSOCIATED(wf_buf2_cptr)) then
    call C_F_POINTER (wf_buf2_cptr, wf_buf2, [nx, ny])
  else
    cache_ok = .false.
  endif
endif

if (.not. cache_ok) then
  call out_io (s_error$, r_name, &
        'FFTW COULD NOT ALLOCATE A WORK BUFFER OR CREATE A PLAN FOR A \i0\ BY \i0\ TRANSFORM.', &
        i_array = [nx, ny])
  return
endif

buf => wf_buf
if (present(buf2)) buf2 => wf_buf2
if (present(err_flag)) err_flag = .false.

end subroutine wavefront_fft2_buffers

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_fft2_aligned (which, direction, err_flag)
!
! Routine to transform, where it lies, one of the calling thread's work buffers
! (wavefront_fft2_buffers). The sign convention and the normalization are those of
! wavefront_fft2.
!
! The buffer is named by its number and is not an argument. A plan run on an array
! whose alignment is not the one it was made with is undefined, and gfortran passes a
! packed copy of a pointer array to a contiguous dummy, so an array argument here would
! be transformed in a temporary.
!
! Input:
!   which      -- integer: 1 for the first buffer, 2 for the second.
!   direction  -- integer: wf_fft_forward$ or wf_fft_backward$.
!
! Output:
!   err_flag   -- logical: Set True if the buffer is not allocated or the direction is
!                   not one of the two.
!-

subroutine wavefront_fft2_aligned (which, direction, err_flag)

include 'fftw3.f03'

integer which, direction
logical err_flag
character(*), parameter :: r_name = 'wavefront_fft2_aligned'

!

err_flag = .true.

if (direction /= wf_fft_forward$ .and. direction /= wf_fft_backward$) then
  call out_io (s_error$, r_name, 'BAD TRANSFORM DIRECTION: \i0\ ', i_array = [direction])
  return
endif

if (wf_cache_nx == 0 .or. (which == 1 .and. .not. associated(wf_buf)) .or. &
    (which == 2 .and. .not. associated(wf_buf2)) .or. which < 1 .or. which > 2) then
  call out_io (s_error$, r_name, 'WORK BUFFER \i0\ OF THIS THREAD IS NOT ALLOCATED.', i_array = [which])
  return
endif

! Passing the same array as both in and out is FFTW's documented way to ask for an
! in-place transform, as bmad/space_charge/fft_interface_mod.f90 does. gfortran -Wall
! warns about it, because the fftw3.f03 interface declares both dummies intent(out).
! Aliasing those would be nonconforming for a normal Fortran procedure. This one is
! bind(C) with assumed-size dummies, so the arguments are passed by address with no copy
! in or out, and the aliasing is safe. Bmad's own build does not enable that warning.
! The second buffer has the first's size and alignment, which is what a plan made on
! the first requires of another array.

if (which == 1) then
  if (direction == wf_fft_forward$) then
    call fftw_execute_dft (wf_plan_fwd, wf_buf, wf_buf)
  else
    call fftw_execute_dft (wf_plan_bwd, wf_buf, wf_buf)
  endif
else
  if (direction == wf_fft_forward$) then
    call fftw_execute_dft (wf_plan_fwd, wf_buf2, wf_buf2)
  else
    call fftw_execute_dft (wf_plan_bwd, wf_buf2, wf_buf2)
  endif
endif

err_flag = .false.

end subroutine wavefront_fft2_aligned

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
integer(C_INT) n_ok, flags
logical cache_ok

!

!$OMP CRITICAL (wavefront_fft_plan_lock)

if (nx /= wf_cache_nx .or. ny /= wf_cache_ny .or. wf_cache_mode /= wf_plan_mode) then
  if (C_ASSOCIATED(wf_plan_fwd)) call fftw_destroy_plan(wf_plan_fwd)
  if (C_ASSOCIATED(wf_plan_bwd)) call fftw_destroy_plan(wf_plan_bwd)
  if (C_ASSOCIATED(wf_buf_cptr)) call fftw_free(wf_buf_cptr)
  if (C_ASSOCIATED(wf_buf2_cptr)) call fftw_free(wf_buf2_cptr)
  wf_buf2_cptr = C_NULL_PTR
  wf_buf2 => null()
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
    ! The plan is one every run and every thread makes alike. FFTW_MEASURE picks the
    ! algorithm by timing candidate transforms, so two runs of the same binary can pick
    ! differently. Every transform after that differs at the ulp level, enough to change a
    ! run's field diagnostics byte for byte. Under wf_plan_wisdom$ the plan is the one
    ! the wisdom text records, a measured plan chosen once and written down. The text is
    ! read, the plan is made from it alone (FFTW_WISDOM_ONLY, which returns null where
    ! the text does not hold this problem), and the text is forgotten again, so no other
    ! plan of the program is made from it. FFTW_PATIENT is the level the text was
    ! measured at. Where it gives no plan, and under wf_plan_estimate$, the plan is
    ! FFTW_ESTIMATE's, a pure function of the problem.
    wf_cache_made = wf_plan_mode

    if (wf_plan_mode == wf_plan_wisdom$) then
      n_ok = fftw_import_wisdom_from_string(wavefront_wisdom_text())
      flags = ior(FFTW_PATIENT, FFTW_WISDOM_ONLY)
      wf_plan_fwd = fftw_plan_dft_2d(ny, nx, wf_buf, wf_buf, FFTW_FORWARD,  flags)
      wf_plan_bwd = fftw_plan_dft_2d(ny, nx, wf_buf, wf_buf, FFTW_BACKWARD, flags)
      call fftw_forget_wisdom()
      if (.not. C_ASSOCIATED(wf_plan_fwd) .or. .not. C_ASSOCIATED(wf_plan_bwd)) then
        if (C_ASSOCIATED(wf_plan_fwd)) call fftw_destroy_plan(wf_plan_fwd)
        if (C_ASSOCIATED(wf_plan_bwd)) call fftw_destroy_plan(wf_plan_bwd)
        wf_cache_made = wf_plan_estimate$
      endif

      ! The text of a file, where the source's text gave no plan: read, planned from
      ! alone and forgotten in the same way.

      if (wf_cache_made == wf_plan_estimate$ .and. allocated(wf_file_wisdom)) then
        n_ok = fftw_import_wisdom_from_string(wf_file_wisdom)
        wf_plan_fwd = fftw_plan_dft_2d(ny, nx, wf_buf, wf_buf, FFTW_FORWARD,  flags)
        wf_plan_bwd = fftw_plan_dft_2d(ny, nx, wf_buf, wf_buf, FFTW_BACKWARD, flags)
        call fftw_forget_wisdom()
        if (C_ASSOCIATED(wf_plan_fwd) .and. C_ASSOCIATED(wf_plan_bwd)) then
          wf_cache_made = wf_plan_wisdom_file$
        else
          if (C_ASSOCIATED(wf_plan_fwd)) call fftw_destroy_plan(wf_plan_fwd)
          if (C_ASSOCIATED(wf_plan_bwd)) call fftw_destroy_plan(wf_plan_bwd)
        endif
      endif
    endif

    if (wf_cache_made /= wf_plan_wisdom$ .and. wf_cache_made /= wf_plan_wisdom_file$) then
      flags = FFTW_ESTIMATE
      if (wf_cache_made == wf_plan_measure$) flags = FFTW_MEASURE
      wf_plan_fwd = fftw_plan_dft_2d(ny, nx, wf_buf, wf_buf, FFTW_FORWARD,  flags)
      wf_plan_bwd = fftw_plan_dft_2d(ny, nx, wf_buf, wf_buf, FFTW_BACKWARD, flags)
    endif

    if (C_ASSOCIATED(wf_plan_fwd) .and. C_ASSOCIATED(wf_plan_bwd)) then
      wf_cache_nx = nx; wf_cache_ny = ny
      wf_cache_mode = wf_plan_mode
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
! Subroutine wavefront_fft2_set_plan (name, err_flag)
!
! Routine to set how the plans of wavefront_fft2 are made: 'wisdom', the plan the text
! of wavefront_wisdom_mod records, else the plan a wisdom file records
! (wavefront_fft2_set_wisdom_file), else the estimated plan. 'wisdom_file' is what
! wavefront_fft2_plan_made answers where the file made the plan, and sets the same,
! 'estimate', FFTW_ESTIMATE, or 'measure', FFTW_MEASURE, whose choice can differ from
! run to run. Call it serially, before any transform. A thread whose cache was filled
! another way plans again at its next transform.
!
! Input:
!   name      -- character(*): 'wisdom', 'wisdom_file', 'estimate' or 'measure'.
!
! Output:
!   err_flag  -- logical: Set True if name is none of the three, False otherwise.
!-

subroutine wavefront_fft2_set_plan (name, err_flag)

character(*) name
logical err_flag
character(*), parameter :: r_name = 'wavefront_fft2_set_plan'

!

err_flag = .false.

select case (name)
case ('wisdom');    wf_plan_mode = wf_plan_wisdom$
case ('wisdom_file');  wf_plan_mode = wf_plan_wisdom$
case ('estimate');  wf_plan_mode = wf_plan_estimate$
case ('measure');   wf_plan_mode = wf_plan_measure$
case default
  call out_io (s_error$, r_name, 'UNKNOWN FFT PLAN: ' // trim(name))
  err_flag = .true.
end select

end subroutine wavefront_fft2_set_plan

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_fft2_plan_made () result (name)
!
! Routine to return how the plans in the calling thread's cache were made: 'wisdom',
! 'wisdom_file', 'estimate' or 'measure'. Under 'wisdom' this is 'wisdom_file' where
! the source's text gave no plan for the cached size and a file's did, and 'estimate'
! where neither did. Every thread plans alike, so after
! wavefront_fft2_plan_threads the calling thread's answer is every thread's.
!
! Output:
!   name  -- character(12): How the plans were made, blank if the cache is empty.
!-

function wavefront_fft2_plan_made () result (name)

character(12) name

!

name = ''
if (wf_cache_nx > 0 .and. wf_cache_made > 0) name = wf_plan_name(wf_cache_made)

end function wavefront_fft2_plan_made

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_fft2_set_wisdom_file (file_name, found)
!
! Routine to read a file of FFTW wisdom, as wavefront_fft2_make_wisdom writes it, for
! the plans of wavefront_fft2 under 'wisdom' (wavefront_fft2_set_plan). The file is
! read here once and never again. A blank name, or a file that is not there or is
! empty, leaves no file's wisdom in use. Call it serially, before any transform.
!
! Input:
!   file_name  -- character(*): The file, or blank for none.
!
! Output:
!   found      -- logical: True if the file was read.
!-

subroutine wavefront_fft2_set_wisdom_file (file_name, found)

character(*) file_name
logical found
integer iu, ios, n

!

found = .false.
if (allocated(wf_file_wisdom)) deallocate (wf_file_wisdom)
if (file_name == '') return

inquire (file = trim(file_name), exist = found, size = n)
if (.not. found .or. n < 1) then
  found = .false.
  return
endif

open (newunit = iu, file = trim(file_name), access = 'stream', form = 'unformatted', action = 'read', &
      status = 'old', iostat = ios)
if (ios /= 0) then
  found = .false.
  return
endif
allocate (character(n+1, kind = C_CHAR) :: wf_file_wisdom)
read (iu, iostat = ios) wf_file_wisdom(1:n)
close (iu)
wf_file_wisdom(n+1:n+1) = C_NULL_CHAR
if (ios /= 0) then
  deallocate (wf_file_wisdom)
  found = .false.
endif

end subroutine wavefront_fft2_set_wisdom_file

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine wavefront_fft2_make_wisdom (file_name, n_side, state, seconds, err_flag)
!
! Routine to measure the plans of wavefront_fft2 for square grids and write FFTW's
! record of them, its wisdom, to a file that wavefront_fft2_set_wisdom_file reads.
! Each grid is planned at FFTW_PATIENT, in place, forward and backward, as
! wavefront_fft2_plan plans it. A grid the text of wavefront_wisdom_mod already gives a
! plan for is not measured, since a run plans it from that text, and a grid the file
! already holds is not measured again. The file is written only if a grid was measured,
! with what it held before, so a second call for the same grids leaves its bytes alone.
!
! A file that is there and that this FFTW does not accept, another version's, is an
! error, and it is not overwritten. Call this serially, outside any transform.
!
! Input:
!   file_name   -- character(*): The file to add to, made if it is not there.
!   n_side(:)   -- integer: The grids' points on a side.
!
! Output:
!   state(:)    -- integer: For each grid, wf_wisdom_in_source$, wf_wisdom_in_file$ or
!                    wf_wisdom_measured$.
!   seconds(:)  -- real(rp): The time the measurement of a grid took, 0 if not measured.
!   err_flag    -- logical: Set True if the file could not be read or written or a plan
!                    could not be made.
!-

subroutine wavefront_fft2_make_wisdom (file_name, n_side, state, seconds, err_flag)

include 'fftw3.f03'

character(*) file_name
integer n_side(:), state(:)
real(rp) seconds(:)
logical err_flag

type(C_PTR) buf_cptr, p_f, p_b
complex(wf_rp), pointer :: buf(:,:)
integer i, n
integer(C_INT) n_ok
integer(8) t_0, t_1, t_rate
logical there, any_new
character(*), parameter :: r_name = 'wavefront_fft2_make_wisdom'

!

err_flag = .true.
state = 0
seconds = 0

! The grids the source's text already gives a plan for.

call fftw_forget_wisdom()
n_ok = fftw_import_wisdom_from_string(wavefront_wisdom_text())
do i = 1, size(n_side)
  if (has_plan(n_side(i))) state(i) = wf_wisdom_in_source$
enddo
call fftw_forget_wisdom()

! The file's wisdom, then each grid it does not hold, measured on top of it.

inquire (file = trim(file_name), exist = there)
if (there) then
  n_ok = fftw_import_wisdom_from_filename(trim(file_name) // C_NULL_CHAR)
  if (n_ok == 0) then
    call out_io (s_error$, r_name, 'THIS FFTW DOES NOT ACCEPT THE WISDOM OF: ' // trim(file_name), &
                 'THE FILE IS LEFT AS IT IS. REMOVE IT, OR NAME ANOTHER.')
    call fftw_forget_wisdom()
    return
  endif
endif

any_new = .false.
call system_clock (t_0, t_rate)
do i = 1, size(n_side)
  if (state(i) /= 0) cycle
  n = n_side(i)
  if (has_plan(n)) then
    state(i) = wf_wisdom_in_file$
    cycle
  endif
  buf_cptr = fftw_alloc_complex(int(n, C_SIZE_T) * int(n, C_SIZE_T))
  if (.not. C_ASSOCIATED(buf_cptr)) exit
  call C_F_POINTER (buf_cptr, buf, [n, n])
  call system_clock (t_0)
  p_f = fftw_plan_dft_2d(n, n, buf, buf, FFTW_FORWARD,  FFTW_PATIENT)
  p_b = fftw_plan_dft_2d(n, n, buf, buf, FFTW_BACKWARD, FFTW_PATIENT)
  call system_clock (t_1)
  if (C_ASSOCIATED(p_f)) call fftw_destroy_plan(p_f)
  if (C_ASSOCIATED(p_b)) call fftw_destroy_plan(p_b)
  call fftw_free(buf_cptr)
  if (.not. C_ASSOCIATED(p_f) .or. .not. C_ASSOCIATED(p_b)) exit
  state(i) = wf_wisdom_measured$
  seconds(i) = real(t_1 - t_0, rp) / t_rate
  any_new = .true.
enddo

if (any(state == 0)) then
  call out_io (s_error$, r_name, 'FFTW COULD NOT MAKE A PLAN FOR A \i0\ BY \i0\ TRANSFORM.', &
               i_array = [n_side(minloc(state, 1)), n_side(minloc(state, 1))])
  call fftw_forget_wisdom()
  return
endif

if (any_new) then
  n_ok = fftw_export_wisdom_to_filename(trim(file_name) // C_NULL_CHAR)
  if (n_ok == 0) then
    call out_io (s_error$, r_name, 'CANNOT WRITE THE WISDOM TO: ' // trim(file_name))
    call fftw_forget_wisdom()
    return
  endif
endif

call fftw_forget_wisdom()
err_flag = .false.

!------------------------------------------------------------------------------
contains

! True if the planner's wisdom holds the in-place plans of an n by n transform.

function has_plan (n) result (held)

integer n
logical held
type(C_PTR) b_cptr, q_f, q_b
complex(wf_rp), pointer :: b(:,:)
integer(C_INT) flags

!

held = .false.
b_cptr = fftw_alloc_complex(int(n, C_SIZE_T) * int(n, C_SIZE_T))
if (.not. C_ASSOCIATED(b_cptr)) return
call C_F_POINTER (b_cptr, b, [n, n])
flags = ior(FFTW_PATIENT, FFTW_WISDOM_ONLY)
q_f = fftw_plan_dft_2d(n, n, b, b, FFTW_FORWARD,  flags)
q_b = fftw_plan_dft_2d(n, n, b, b, FFTW_BACKWARD, flags)
held = C_ASSOCIATED(q_f) .and. C_ASSOCIATED(q_b)
if (C_ASSOCIATED(q_f)) call fftw_destroy_plan(q_f)
if (C_ASSOCIATED(q_b)) call fftw_destroy_plan(q_b)
call fftw_free(b_cptr)

end function has_plan

end subroutine wavefront_fft2_make_wisdom

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_fftw_version () result (version)
!
! Routine to return the version FFTW states in the first word of its wisdom, as
! fftw-3.3.11, which names the wisdom a file can hold.
!
! Output:
!   version  -- character(32): The version, blank if FFTW exported nothing.
!-

function wavefront_fftw_version () result (version)

include 'fftw3.f03'

character(32) version
type(C_PTR) text_cptr
character(kind = C_CHAR), pointer :: text(:)
integer i

interface
  subroutine c_free (ptr) bind(C, name = 'free')
    import :: C_PTR
    type(C_PTR), value :: ptr
  end subroutine
end interface

!

version = ''
text_cptr = fftw_export_wisdom_to_string()
if (.not. C_ASSOCIATED(text_cptr)) return
call C_F_POINTER (text_cptr, text, [40])

! The text opens with a parenthesis, and the version runs to the first blank.

do i = 2, 33
  if (text(i) == ' ' .or. text(i) == C_NULL_CHAR) exit
  version(i-1:i-1) = text(i)
enddo
call c_free (text_cptr)

end function wavefront_fftw_version

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function wavefront_wisdom_text () result (text)
!
! Routine to return the wisdom text of wavefront_wisdom_mod as one C string, its lines
! joined and a null at the end, as fftw_import_wisdom_from_string takes it.
!
! Output:
!   text  -- character(:), allocatable: The wisdom text.
!-

function wavefront_wisdom_text () result (text)

character(:, kind = C_CHAR), allocatable :: text
integer i

!

text = ''
do i = 1, size(wf_wisdom_line)
  text = text // trim(wf_wisdom_line(i)) // new_line('a')
enddo
text = text // C_NULL_CHAR

end function wavefront_wisdom_text

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
if (C_ASSOCIATED(wf_buf2_cptr)) call fftw_free(wf_buf2_cptr)
wf_buf2_cptr = C_NULL_PTR
wf_buf2 => null()
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
! own. It is the quantity docs/startup-noise.md splits on. Inside it lies the radiation
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
