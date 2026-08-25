!+
! Module fel_beam_mod
!
! This module holds the packed particle representation for FEL tracking, its Genesis 1.3
! Version 4 particle dump I/O, its conversion to and from Bmad's coord_struct at element
! boundaries, and the per-slice beam diagnostics.
!
! Coordinates: Bmad's, exactly. Each slice carries structure-of-arrays copies of
! coord_struct%vec(1:6) plus a weight:
!
!   x, y      [m]
!   px, py    transverse momentum over the reference momentum, P_xy/p0
!   z         -beta*c*(t - t_ref)  [m], the Bmad longitudinal coordinate
!   pz        (p - p0)/p0
!   weight    macroparticle charge [C], mapping to coord_struct%charge
!
! so conversion at the seam is a plain copy. (z, pz) is the conjugate pair of Bmad's
! s-based Hamiltonian. The normalization reference is stored once, as p0c [eV] -- Bmad's
! own convention, asserted against each element's p0c at conversion time. The equivalent
! quantities Genesis works in are derived, never stored: p0/(m_e c) from fel_p0_mc and
! the reference gamma from fel_gamma0. (Genesis carries its reference gamma, gammaref, as
! an independent run parameter, and importing a dump defines p0c from it. Deriving gamma
! back from p0c reproduces gammaref to ~1 ulp, far below the 8e-7 constants floor of any
! Genesis comparison.)
!
! The ponderomotive phase is derived, not stored. Genesis's per-particle theta is the sum
! of a common reference advance (the undulator's ku term, and the drift slippage term
! ks/(2 gamma0^2)) and a particle-specific lag. The common part lives in one scalar per
! beam, phi0, maintained by the tracker. The particle part is the Bmad z:
!
!   theta_j = phi0 - ks * tau_j,     tau_j = -z_j / beta_j = c*(t_j - t_ref)
!
! This reference-offset split is the FP32-safe formulation, and it removes the theta-wrap
! hazard outright: z does not wrap, so there is no theta-wrap-plus-slice-index update to
! get wrong at slice migration. The slice index is derived from z when needed. (A future
! single precision GPU struct still needs per-slice re-referencing of z, since the phase
! needs ~1e-6 rad across ~1e5 wavelengths of bunch. That is a choice deferred to the
! device struct.)
!
! Weights are carried from day one: every reduction here and in
! fel_track_mod is weighted. Slice current is derived as I = c * sum(w) / slice_spacing, and
! N_eff = (sum w)^2 / sum w^2 is a per-slice diagnostic. A uniform-weight beam reproduces
! Genesis, which the benchmark checks.
!
! Why packed arrays at all: the FEL step advances every particle every internal step, and
! coord_struct is ~224 bytes against the ~56 needed. coord_struct appears only at element
! boundaries. The arrays are allocated to a capacity that may exceed the fill count n.
! Slice migration can then move particles without per-step reallocation (fel-physics.tex
! sec:migration).
!-

module fel_beam_mod

use bmad
use hdf5_interface

implicit none

! All physical constants come from sim_utils (m_electron, c_light, mu_0_vac). Genesis
! carries its own values (EG an impedance of free space truncated to 376.73 where
! mu_0_vac*c_light is 376.7303...). During the original transcription validation this
! module transcribed them to get transcription-level agreement. That validation is banked
! (lucifer/README.md). The code now uses Bmad's constants, and the ~8e-7 relative
! difference against Genesis is the accepted comparison floor.

!+
! Struct fel_slice_struct
!
! One beam slice in packed structure-of-arrays form, Bmad coordinates plus weight.
! Live particles are 1:n. The arrays may be larger.
!-

type fel_slice_struct
  real(rp), allocatable :: x(:)        ! [m]
  real(rp), allocatable :: px(:)       ! P_x/p0
  real(rp), allocatable :: y(:)        ! [m]
  real(rp), allocatable :: py(:)       ! P_y/p0
  real(rp), allocatable :: z(:)        ! -beta*c*(t - t_ref) [m]
  real(rp), allocatable :: pz(:)       ! (p - p0)/p0
  real(rp), allocatable :: weight(:)   ! Macroparticle charge [C]
  integer :: n = 0                     ! Fill count.
end type

!+
! Struct fel_beam_struct
!
! The whole beam: slices, the normalization reference, the common phase, and the dump
! metadata.
!-

type fel_beam_struct
  type (fel_slice_struct), allocatable :: slice(:)
  real(rp) :: p0c = 0              ! Reference momentum times c [eV]. The single stored
                                   !   reference. Derive p0/(m_e c) and the reference gamma
                                   !   with fel_p0_mc and fel_gamma0.
  real(rp) :: phi0 = 0             ! Common ponderomotive reference phase [rad].
  real(rp) :: wavelength = 0       ! Radiation wavelength [m]. 'slicelength' in a Genesis dump.
  real(rp) :: slice_spacing = 0    ! Longitudinal slice spacing [m]. 'slicespacing' in a Genesis dump.
  real(rp) :: s0 = 0               ! Start of the time window [m]. 'refposition' in a Genesis dump.
  integer :: nbins = 0             ! Beamlet size at generation. Carried for dump round trips.
  logical :: one4one = .false.     ! Genesis one4one flag. Carried for dump round trips.
  logical :: quiver_in_px = .false. ! Momentum convention flag (fel-physics.tex
                                   !   sec:unaveraged): False = stored px excludes the
                                   !   undulator quiver (the averaged chart, and the only
                                   !   convention any averaged/seam code accepts). True =
                                   !   the unaveraged tracker is mid-region and px is the
                                   !   full kinetic momentum. Asserted at region edges: a
                                   !   hard-edge handoff injects K/gamma of spurious px.
end type

!+
! Struct fel_slice_diag_struct
!
! Per-slice beam diagnostics at one output position: the quantities of Genesis's
! DiagBeam::calc, weighted, under Bmad-style names and normalizations. mean_/sigma_
! follow the convention of wavefront_transverse_moments. gamma is the Lorentz factor,
! named as such because Bmad's 'energy' is the total energy in eV, which none of these
! are. The Genesis output dataset each one compares against is noted.
!
! Scope, by decision: this struct is the Genesis-comparison instrument and stays limited
! to quantities a Genesis output can check. Production moment diagnostics (mean vector,
! sigma matrix, emittances) go through Bmad's bunch_params_struct at the seam, where
! calc_bunch_params already owns the definitions.
!-

type fel_slice_diag_struct
  real(rp) :: mean_energy = 0      ! Weighted <E> = <gamma>*m_e c^2 [eV], Bmad's
                                   !   convention (energy means eV). Genesis 'energy'
                                   !   is <gamma>. Divide by m_electron to compare.
  real(rp) :: sigma_energy = 0     ! rms energy spread [eV]. Genesis 'energyspread'
                                   !   is the same in units of gamma.
  real(rp) :: bunching = 0         ! |sum w e^{i theta}| / sum w. Genesis 'bunching'.
  real(rp) :: bunching_phase = 0   ! arg(sum w e^{i theta}). Genesis 'bunchingphase'.
  real(rp) :: mean_x = 0, mean_y = 0    ! Weighted centroid [m]. Genesis 'xposition', 'yposition'.
  real(rp) :: sigma_x = 0, sigma_y = 0  ! Weighted rms size [m]. Genesis 'xsize', 'ysize'.
  real(rp) :: mean_px = 0, mean_py = 0  ! Weighted <P_x/p0>, <P_y/p0>, Bmad normalization.
                                        !   Genesis 'pxposition' is this times p0/(m_e c).
  real(rp) :: n_eff = 0            ! (sum w)^2 / sum w^2 (no Genesis counterpart).
  real(rp) :: current = 0          ! c * sum(w) / slice_spacing [A]. Genesis 'current'.
end type

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_gamma_of (p0_mc, pz) result (gamma)
!
! Routine to return a particle's Lorentz factor, derived from the stored momentum
! coordinates. Exact algebra: gamma = sqrt(P^2 + 1) with P = p0_mc * (1 + pz) the total
! momentum over m_e c.
!
! Input:
!   p0_mc   -- real(rp): Reference momentum p0/(m_e c). See fel_p0_mc.
!   pz      -- real(rp): Bmad pz = (p - p0)/p0.
!
! Output:
!   gamma   -- real(rp): Lorentz factor.
!-

elemental function fel_gamma_of (p0_mc, pz) result (gamma)
real(rp), intent(in) :: p0_mc, pz
real(rp) gamma
gamma = sqrt((p0_mc * (1 + pz))**2 + 1)
end function

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_beta_of (p0_mc, pz) result (beta)
!
! Routine to return a particle's velocity beta, derived from the stored momentum
! coordinates. Exact algebra: beta = P/gamma = P/sqrt(P^2 + 1) with P = p0_mc * (1 + pz)
! the total momentum over m_e c. The trajectory time follows as tau = -z/beta, which is
! c*(t - t_ref) in meters.
!
! Input:
!   p0_mc   -- real(rp): Reference momentum p0/(m_e c). See fel_p0_mc.
!   pz      -- real(rp): Bmad pz = (p - p0)/p0.
!
! Output:
!   beta    -- real(rp): Velocity over c.
!-

elemental function fel_beta_of (p0_mc, pz) result (beta)
real(rp), intent(in) :: p0_mc, pz
real(rp) beta, p_mc
p_mc = p0_mc * (1 + pz)
beta = p_mc / sqrt(p_mc**2 + 1)
end function

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_p0_mc (beam) result (p0_mc)
!
! Routine to return the dimensionless reference momentum p0/(m_e c) = p0c/m_electron,
! derived from the beam's one stored reference p0c. Bmad px times this is Genesis px.
! Hoist this out of particle loops. It is one division, but there is no reason to pay
! it per particle.
!
! Input:
!   beam    -- fel_beam_struct: The beam.
!
! Output:
!   p0_mc   -- real(rp): Reference momentum over m_e c.
!-

elemental function fel_p0_mc (beam) result (p0_mc)
type (fel_beam_struct), intent(in) :: beam
real(rp) p0_mc
p0_mc = beam%p0c / m_electron
end function

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_gamma0 (beam) result (gamma0)
!
! Routine to return the gamma of the reference momentum (Genesis's gammaref), derived
! from the beam's one stored reference p0c. Hoist this out of particle loops. It is one
! sqrt, but there is no reason to pay it per particle.
!
! Input:
!   beam    -- fel_beam_struct: The beam.
!
! Output:
!   gamma0  -- real(rp): Reference Lorentz factor sqrt((p0c/m_electron)^2 + 1).
!-

elemental function fel_gamma0 (beam) result (gamma0)
type (fel_beam_struct), intent(in) :: beam
real(rp) gamma0
gamma0 = sqrt((beam%p0c / m_electron)**2 + 1)
end function

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_assert_averaged_chart (beam, where_str, err_flag)
!
! Routine to refuse by name when the beam's stored px still carries the undulator quiver
! (the unaveraged tracker left its region without restoring the averaged convention).
! Every averaged-physics or seam entry point calls this: a hard-edge handoff injects
! spurious transverse momentum of order K/gamma -- comparable to the beam divergence --
! and presents downstream as unexplained emittance growth (fel-physics.tex
! sec:unaveraged).
!
! Input:
!   beam        -- fel_beam_struct: The beam whose chart flag is checked.
!   where_str   -- character(*): Caller identification for the error message.
!
! Output:
!   err_flag    -- logical: Set True if the beam is in the quiver chart, False otherwise.
!-

subroutine fel_assert_averaged_chart (beam, where_str, err_flag)

type (fel_beam_struct) beam
character(*) where_str
logical err_flag
character(*), parameter :: r_name = 'fel_assert_averaged_chart'

!

err_flag = .false.
if (.not. beam%quiver_in_px) return

err_flag = .true.
call out_io (s_error$, r_name, 'BEAM px STILL CARRIES THE UNDULATOR QUIVER (quiver_in_px IS SET):', &
        'THE UNAVERAGED TRACKER DID NOT HAND BACK THE AVERAGED CONVENTION AT A RAMP END.', &
        'AVERAGED PHYSICS AND THE SEAM REFUSE THIS CHART. AT: ' // trim(where_str))

end subroutine fel_assert_averaged_chart

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_theta (beam, sl, ip, ks) result (theta)
!
! Routine to return the ponderomotive phase of particle ip: theta = phi0 - ks*tau with
! tau = -z/beta. See the module header for the reference-offset split.
!
! Input:
!   beam    -- fel_beam_struct: The beam (for phi0 and the reference momentum).
!   sl      -- fel_slice_struct: Slice holding the particle.
!   ip      -- integer: Particle index within the slice.
!   ks      -- real(rp): Radiation wavenumber twopi/wavelength [1/m].
!
! Output:
!   theta   -- real(rp): Ponderomotive phase [rad].
!-

function fel_theta (beam, sl, ip, ks) result (theta)

type (fel_beam_struct) beam
type (fel_slice_struct) sl
integer ip
real(rp) ks, theta, beta

!

beta = fel_beta_of(fel_p0_mc(beam), sl%pz(ip))
theta = beam%phi0 + ks * sl%z(ip) / beta      ! phi0 - ks*(-z/beta)

end function fel_theta

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_read_genesis4_beam (beam, file_name, gamma0, err_flag)
!
! Routine to read a Genesis 1.3 Version 4 particle dump (.par.h5) into a packed beam,
! converting Genesis coordinates (x, y in m; px, py = gamma*beta; theta; gamma) to Bmad's.
!
! The conversion, exact:
!
!   px_bmad = px_genesis / p0_mc                    (p0_mc = sqrt(gamma0^2 - 1))
!   pz      = (sqrt(gamma^2 - 1) - p0_mc) / p0_mc
!   tau     = -theta / ks,   z = -beta * tau        (phi0 starts at zero)
!   weight  = I * slice_spacing / (c * n)  [C]      (uniform: the dump carries no weights)
!
! ks comes from the dump's own wavelength. p0c is set from gamma0 and Bmad's electron
! mass, which is what makes the Bmad-side tracking see exactly this normalization.
!
! Input:
!   file_name   -- character(*): File to read.
!   gamma0      -- real(rp): Genesis's reference gamma for the run.
!
! Output:
!   beam        -- fel_beam_struct: Beam read from the file.
!   err_flag    -- logical: Set True on error, False otherwise.
!-

subroutine fel_read_genesis4_beam (beam, file_name, gamma0, err_flag)

type (fel_beam_struct), target :: beam
type (fel_slice_struct), pointer :: sl
integer(hid_t) f_id, g_id
integer is, ip, n_slice, np, ivec(1), h5_err
real(rp) rvec(1), gamma0, ks, current, w_uniform, gam, p_mc, beta, tau, p0_mc
real(rp), allocatable :: work_gamma(:), work_theta(:)
logical err_flag, err
character(*) file_name
character(*), parameter :: r_name = 'fel_read_genesis4_beam'
character(20) group_name
type (hdf5_info_struct) info

!

err_flag = .true.

if (gamma0 <= 1) then
  call out_io (s_error$, r_name, 'gamma0 MUST EXCEED 1. GOT: \es12.4\ ', r_array = [gamma0])
  return
endif

call hdf5_open_file (file_name, 'READ', f_id, err);  if (err) return

call hdf5_read_dataset_int (f_id, 'slicecount', ivec, err, 'slicecount');  if (err) return
n_slice = ivec(1)
if (n_slice < 1) then
  call out_io (s_error$, r_name, 'FILE HAS A NON-POSITIVE slicecount: \i0\ ', i_array = [n_slice])
  return
endif

call hdf5_read_dataset_int (f_id, 'beamletsize', ivec, err, 'beamletsize');  if (err) return
beam%nbins = ivec(1)

call hdf5_read_dataset_int (f_id, 'one4one', ivec, err, 'one4one');  if (err) return
beam%one4one = (ivec(1) /= 0)

call hdf5_read_dataset_real (f_id, 'slicelength', rvec, err, 'slicelength');  if (err) return
beam%wavelength = rvec(1)

call hdf5_read_dataset_real (f_id, 'slicespacing', rvec, err, 'slicespacing');  if (err) return
beam%slice_spacing = rvec(1)

call hdf5_read_dataset_real (f_id, 'refposition', rvec, err, 'refposition');  if (err) return
beam%s0 = rvec(1)

! p0c is the one stored reference: the momentum whose gamma is Genesis's gammaref.

p0_mc = sqrt(gamma0**2 - 1)
beam%p0c = p0_mc * m_electron
beam%phi0 = 0
ks = twopi / beam%wavelength

if (allocated(beam%slice)) deallocate(beam%slice)
allocate (beam%slice(n_slice))

do is = 1, n_slice
  sl => beam%slice(is)
  write (group_name, '(a, i0.6)') 'slice', is
  g_id = hdf5_open_group (f_id, trim(group_name), err, .true.);  if (err) return

  ! Particle count from the actual dataset extent: H5LT reads have no buffer bound, so
  ! extents are checked, never assumed.

  info = hdf5_object_info (g_id, 'gamma', err, .true.);  if (err) return
  np = int(info%data_dim(1))

  call fel_slice_reallocate (sl, np)
  sl%n = np

  call hdf5_read_dataset_real (g_id, 'current', rvec, err, trim(group_name) // '/current');  if (err) return
  current = rvec(1)

  allocate (work_gamma(np), work_theta(np))
  call hdf5_read_dataset_real (g_id, 'gamma', work_gamma, err, trim(group_name) // '/gamma');  if (err) return
  call hdf5_read_dataset_real (g_id, 'theta', work_theta, err, trim(group_name) // '/theta');  if (err) return
  call hdf5_read_dataset_real (g_id, 'x',  sl%x(1:np),  err, trim(group_name) // '/x');   if (err) return
  call hdf5_read_dataset_real (g_id, 'y',  sl%y(1:np),  err, trim(group_name) // '/y');   if (err) return
  call hdf5_read_dataset_real (g_id, 'px', sl%px(1:np), err, trim(group_name) // '/px');  if (err) return
  call hdf5_read_dataset_real (g_id, 'py', sl%py(1:np), err, trim(group_name) // '/py');  if (err) return

  ! Convert in place: px, py from gamma*beta to P/p0, and (theta, gamma) to (z, pz).

  w_uniform = 0
  if (np > 0) w_uniform = current * beam%slice_spacing / (c_light * np)

  do ip = 1, np
    gam = work_gamma(ip)
    p_mc = sqrt(gam**2 - 1)
    beta = p_mc / gam

    sl%px(ip) = sl%px(ip) / p0_mc
    sl%py(ip) = sl%py(ip) / p0_mc
    sl%pz(ip) = (p_mc - p0_mc) / p0_mc

    tau = -work_theta(ip) / ks              ! theta = phi0 - ks*tau with phi0 = 0.
    sl%z(ip) = -beta * tau

    sl%weight(ip) = w_uniform
  enddo
  deallocate (work_gamma, work_theta)

  call H5Gclose_f (g_id, h5_err)
enddo

call H5Fclose_f (f_id, h5_err)

err_flag = .false.

end subroutine fel_read_genesis4_beam

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_read_openpmd_beam (beam, file_name, gamma0, ele, err_flag)
!
! Routine to read an openPMD-beamphysics particle file written by
! fel_write_openpmd_beam back into a packed beam, weights included. The inverse of that
! routine: one bunch per nonempty slice, placed by the file's felSliceIndex, with the
! window's empty slices restored empty.
!
! A file without the fel* window attributes is refused by name. Such a file describes a
! bunch, not an FEL time window, and belongs on the import path (dist_file), which
! resamples it into slices instead of assuming a slicing it does not carry.
!
! Input:
!   file_name   -- character(*): File to read.
!   gamma0      -- real(rp): Reference gamma for the run, from the lattice.
!   ele         -- ele_struct: Element to convert the coordinates against.
!
! Output:
!   beam        -- fel_beam_struct: Beam read from the file.
!   err_flag    -- logical: Set True on error, False otherwise.
!-

subroutine fel_read_openpmd_beam (beam, file_name, gamma0, ele, err_flag)

type (fel_beam_struct), target :: beam
type (fel_slice_struct), pointer :: sl
type (ele_struct) ele
type (beam_struct) beam_b
integer(hid_t) f_id
integer ib, is, np, nb, n_slice, h5_err, ivec(1)
integer, allocatable :: slice_ix(:)
real(rp) gamma0, rvec(1)
logical err_flag, err
character(*) file_name
character(*), parameter :: r_name = 'fel_read_openpmd_beam'

!

err_flag = .true.

if (gamma0 <= 1) then
  call out_io (s_error$, r_name, 'gamma0 MUST EXCEED 1. GOT: \es12.4\ ', r_array = [gamma0])
  return
endif

! The window attributes first: without them the file cannot be placed into slices.

call hdf5_open_file (file_name, 'READ', f_id, err);  if (err) return

n_slice = -1
call hdf5_read_attribute_int (f_id, 'felNumSlice', n_slice, err, .false.)
if (err .or. n_slice == -1) then
  call out_io (s_error$, r_name, &
      'openPMD PARTICLE FILE CARRIES NO FEL WINDOW ATTRIBUTES: ' // trim(file_name), &
      'IT DESCRIBES A BUNCH, NOT A SLICED TIME WINDOW, SO ITS SLICING IS UNKNOWN.', &
      'READ IT WITH dist_file, WHICH RESAMPLES A BUNCH INTO SLICES.')
  return
endif

if (n_slice < 1) then
  call out_io (s_error$, r_name, 'FILE HAS A NON-POSITIVE felNumSlice: \i0\ ', i_array = [n_slice])
  return
endif

call re_allocate (slice_ix, n_slice)
call hdf5_read_attribute_int (f_id, 'felSliceIndex', slice_ix, err, .true.);  if (err) return
call hdf5_read_attribute_int (f_id, 'felNbins', beam%nbins, err, .true.);  if (err) return
call hdf5_read_attribute_int (f_id, 'felOne4one', ivec(1), err, .true.);  if (err) return
beam%one4one = (ivec(1) /= 0)
call hdf5_read_attribute_real (f_id, 'felSliceSpacing', beam%slice_spacing, err, .true.);  if (err) return
call hdf5_read_attribute_real (f_id, 'felWavelength', beam%wavelength, err, .true.);  if (err) return
call hdf5_read_attribute_real (f_id, 'felWindowStart', beam%s0, err, .true.);  if (err) return
call H5Fclose_f (f_id, h5_err)

! p0c is the one stored reference, as on the Genesis path. phi0 restarts at zero: the
! particle lag lives in z, so the reference phase is free.

beam%p0c = sqrt(gamma0**2 - 1) * m_electron
beam%phi0 = 0

call hdf5_read_beam (file_name, beam_b, err, ele)
if (err) return

nb = size(beam_b%bunch)
if (nb > n_slice) then
  call out_io (s_error$, r_name, 'FILE HAS MORE BUNCHES THAN SLICES: \i0\ VS \i0\ ', &
               'FILE: ' // trim(file_name), i_array = [nb, n_slice])
  return
endif

if (allocated(beam%slice)) deallocate(beam%slice)
allocate (beam%slice(n_slice))

do is = 1, n_slice
  call fel_slice_reallocate (beam%slice(is), 0)
  beam%slice(is)%n = 0
enddo

do ib = 1, nb
  is = slice_ix(ib)
  if (is < 1 .or. is > n_slice) then
    call out_io (s_error$, r_name, 'FILE felSliceIndex \i0\ IS OUTSIDE THE WINDOW OF \i0\ SLICES.', &
                 i_array = [is, n_slice])
    return
  endif
  sl => beam%slice(is)
  np = size(beam_b%bunch(ib)%particle)
  call fel_slice_reallocate (sl, np)
  sl%n = np
  call fel_bunch_to_slice (beam_b%bunch(ib), ele, sl, err)
  if (err) return
enddo

err_flag = .false.

end subroutine fel_read_openpmd_beam

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_write_openpmd_beam (beam, ele, file_name, err_flag)
!
! Routine to write a packed beam as an openPMD-beamphysics particle file, one bunch per
! slice through Bmad's own hdf5_write_beam. Unlike a Genesis dump this carries the
! PER-PARTICLE weight, as openPMD's macro-charge, so a weighted beam survives the round
! trip.
!
! An empty slice is skipped rather than written: hdf5_write_beam takes the species from
! p(1), which a zero-particle bunch does not have. The window geometry therefore rides
! in root attributes instead of being implied by the file's bunch count, and
! felSliceIndex says which slice each written bunch came from. Reading them back
! restores the empty slices exactly.
!
! Input:
!   beam        -- fel_beam_struct: Beam to write.
!   ele         -- ele_struct: Element the beam sits at, for the coord_struct conversion.
!   file_name   -- character(*): File to create.
!
! Output:
!   err_flag    -- logical: Set True on error, False otherwise.
!-

subroutine fel_write_openpmd_beam (beam, ele, file_name, err_flag)

type (fel_beam_struct), target :: beam
type (fel_slice_struct), pointer :: sl
type (ele_struct) ele
type (beam_struct) beam_b
integer(hid_t) f_id
integer is, nb, nslice, h5_err
integer, allocatable :: slice_ix(:)
logical err_flag, err
character(*) file_name
character(*), parameter :: r_name = 'fel_write_openpmd_beam'

!

err_flag = .true.
nslice = size(beam%slice)

nb = count(beam%slice(1:nslice)%n > 0)
if (nb == 0) then
  call out_io (s_error$, r_name, 'BEAM HAS NO PARTICLES IN ANY SLICE; NOTHING TO WRITE.', &
               'FILE: ' // trim(file_name))
  return
endif

call reallocate_beam (beam_b, nb)
call re_allocate (slice_ix, nb)

nb = 0
do is = 1, nslice
  sl => beam%slice(is)
  if (sl%n == 0) cycle
  nb = nb + 1
  slice_ix(nb) = is
  call fel_slice_to_bunch (beam, sl, ele, beam_b%bunch(nb), err)
  if (err) return
enddo

call hdf5_write_beam (file_name, beam_b%bunch(1:nb), .false., err)
if (err) return

! The window the bunches came from. openPMD readers ignore attributes they do not know,
! and Bmad's reader ignores everything outside its base path.

call hdf5_open_file (file_name, 'APPEND', f_id, err);  if (err) return
call hdf5_write_attribute_int (f_id, 'felSliceIndex', slice_ix(1:nb), err);       if (err) return
call hdf5_write_attribute_int (f_id, 'felNumSlice', nslice, err);                 if (err) return
call hdf5_write_attribute_int (f_id, 'felNbins', beam%nbins, err);                if (err) return
call hdf5_write_attribute_real (f_id, 'felSliceSpacing', beam%slice_spacing, err); if (err) return
call hdf5_write_attribute_real (f_id, 'felWavelength', beam%wavelength, err);      if (err) return
call hdf5_write_attribute_real (f_id, 'felWindowStart', beam%s0, err);             if (err) return
call hdf5_write_attribute_int (f_id, 'felOne4one', merge(1, 0, beam%one4one), err); if (err) return
call H5Fclose_f (f_id, h5_err)

err_flag = .false.

end subroutine fel_write_openpmd_beam

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_write_genesis4_beam (beam, file_name, err_flag)
!
! Routine to write a packed beam as a Genesis 1.3 Version 4 particle dump, converting
! back from Bmad coordinates. Inverse of fel_read_genesis4_beam: theta = phi0 - ks*tau
! reconstructs Genesis's unwrapped theta including the accumulated common phase, gamma
! from pz, px py rescaled by p0_mc, current from the weights. The dump format carries no
! per-particle weight, so a nonuniform-weight beam is REFUSED BY NAME here rather than
! written and silently read back uniform. That limit is the format's, not this
! representation's: fel_write_openpmd_beam carries the weights.
!
! Input:
!   beam        -- fel_beam_struct: Beam to write.
!   file_name   -- character(*): File to create.
!
! Output:
!   err_flag    -- logical: Set True on error, False otherwise.
!-

subroutine fel_write_genesis4_beam (beam, file_name, err_flag)

type (fel_beam_struct), target :: beam
type (fel_slice_struct), pointer :: sl
integer(hid_t) f_id, g_id
integer is, ip, h5_err, one4one_int
real(rp) ks, gam, beta, p_mc, p0_mc, wmin, wmax
real(rp), allocatable :: work(:)
logical err_flag, err
character(*) file_name
character(*), parameter :: r_name = 'fel_write_genesis4_beam'
character(20) group_name

!

err_flag = .true.

if (.not. allocated(beam%slice)) then
  call out_io (s_error$, r_name, 'BEAM HAS NO SLICES.')
  return
endif

! This format has one current per slice and no per-particle weight, so a weighted beam
! would come back uniform. Refuse rather than write a beam nothing can read back.

do is = 1, size(beam%slice)
  sl => beam%slice(is)
  if (sl%n < 2) cycle
  wmin = minval(sl%weight(1:sl%n));  wmax = maxval(sl%weight(1:sl%n))
  if (wmax - wmin <= 1e-12_rp * wmax) cycle
  call out_io (s_error$, r_name, &
      'GENESIS FORMAT CANNOT CARRY PER-PARTICLE WEIGHTS, AND THIS BEAM HAS NONUNIFORM WEIGHTS.', &
      'SLICE \i0\ SPREADS \es12.4\ TO \es12.4\ C OVER \i0\ PARTICLES, TOTAL \es12.4\ C.', &
      'THE FORMAT STORES ONE CURRENT PER SLICE, SO A READ-BACK WOULD RETURN A UNIFORM BEAM.', &
      'WRITE openpmd INSTEAD (beam_formats), WHICH CARRIES THE WEIGHTS.', &
      i_array = [is, sl%n], r_array = [wmin, wmax, sum(sl%weight(1:sl%n))])
  return
enddo

ks = twopi / beam%wavelength
p0_mc = fel_p0_mc(beam)

call hdf5_open_file (file_name, 'WRITE', f_id, err);  if (err) return

one4one_int = 0
if (beam%one4one) one4one_int = 1

call hdf5_write_dataset_real (f_id, 'slicelength',  [beam%wavelength],    err);  if (err) return
call hdf5_write_dataset_real (f_id, 'slicespacing', [beam%slice_spacing], err);  if (err) return
call hdf5_write_dataset_real (f_id, 'refposition',  [beam%s0],            err);  if (err) return
call hdf5_write_dataset_int  (f_id, 'beamletsize',  [beam%nbins],       err);  if (err) return
call hdf5_write_dataset_int  (f_id, 'slicecount',   [size(beam%slice)], err);  if (err) return
call hdf5_write_dataset_int  (f_id, 'one4one',      [one4one_int],      err);  if (err) return

do is = 1, size(beam%slice)
  sl => beam%slice(is)
  write (group_name, '(a, i0.6)') 'slice', is
  call H5Gcreate_f (f_id, trim(group_name), g_id, h5_err)
  if (h5_err < 0) then
    call out_io (s_error$, r_name, 'CANNOT CREATE GROUP: ' // trim(group_name))
    return
  endif

  call hdf5_write_dataset_real (g_id, 'current', [c_light * sum(sl%weight(1:sl%n)) / beam%slice_spacing], err)
  if (err) return

  allocate (work(sl%n))

  do ip = 1, sl%n                                             ! gamma
    work(ip) = fel_gamma_of(p0_mc, sl%pz(ip))
  enddo
  call hdf5_write_dataset_real (g_id, 'gamma', work, err);  if (err) return

  do ip = 1, sl%n                                             ! theta = phi0 + ks*z/beta
    work(ip) = beam%phi0 + ks * sl%z(ip) / fel_beta_of(p0_mc, sl%pz(ip))
  enddo
  call hdf5_write_dataset_real (g_id, 'theta', work, err);  if (err) return

  call hdf5_write_dataset_real (g_id, 'x', sl%x(1:sl%n), err);  if (err) return
  call hdf5_write_dataset_real (g_id, 'y', sl%y(1:sl%n), err);  if (err) return

  work = sl%px(1:sl%n) * p0_mc                                ! Back to gamma*beta_x.
  call hdf5_write_dataset_real (g_id, 'px', work, err);  if (err) return
  work = sl%py(1:sl%n) * p0_mc
  call hdf5_write_dataset_real (g_id, 'py', work, err);  if (err) return

  deallocate (work)
  call H5Gclose_f (g_id, h5_err)
enddo

call H5Fclose_f (f_id, h5_err)

err_flag = .false.

end subroutine fel_write_genesis4_beam

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_slice_reallocate (sl, capacity)
!
! Routine to allocate a slice's arrays to at least the given capacity. Existing live
! particles 1:n are preserved. The fill count n is not changed.
!
! Input:
!   sl        -- fel_slice_struct: Slice to size.
!   capacity  -- integer: Minimum array capacity wanted.
!
! Output:
!   sl        -- fel_slice_struct: Slice with arrays of at least the given size.
!-

subroutine fel_slice_reallocate (sl, capacity)

type (fel_slice_struct) sl
integer capacity

!

if (allocated(sl%x)) then
  if (size(sl%x) >= capacity) return
  call re_allocate (sl%x, capacity)
  call re_allocate (sl%px, capacity)
  call re_allocate (sl%y, capacity)
  call re_allocate (sl%py, capacity)
  call re_allocate (sl%z, capacity)
  call re_allocate (sl%pz, capacity)
  call re_allocate (sl%weight, capacity)
else
  allocate (sl%x(capacity), sl%px(capacity), sl%y(capacity), sl%py(capacity), &
            sl%z(capacity), sl%pz(capacity), sl%weight(capacity))
endif

end subroutine fel_slice_reallocate

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_fawley_noise (theta, weight, n, nbins, n_clamp)
!
! Routine to impose Fawley-style shot noise on a quiet-loaded slice, transcribed from
! Genesis's ShotNoise::applyShotNoise and generalized to per-particle weights. Each
! beamlet's amplitude draws on its REAL electron count, the beamlet's charge over e
! (Genesis's slice-uniform ne/mpart for uniform weights). This makes <|b(h)|^2> =
! 1/N_lambda exact for any cross-beamlet weight distribution (fel-physics.tex
! sec:noise). Weights must be uniform WITHIN a beamlet (the quiet cancellation is per
! beamlet). The first particle's weight speaks for its beamlet. Genesis's silent
! nbl < 1 clamp is kept but counted into n_clamp for the caller to report. Kicks
! accumulate from the unperturbed phases, exactly as Genesis's work array does: two
! ran_uniform draws per (harmonic, beamlet), in Genesis's loop order. The routine is
! shared by the built-in loader and the distribution import, so the two stay one
! implementation.
!
! Input:
!   theta(:)    -- real(rp): Quiet-loaded ponderomotive phases [rad], beamlet-contiguous.
!   weight(:)   -- real(rp): Macroparticle charges [C], uniform within each beamlet.
!   n           -- integer: Number of particles in the slice (a multiple of nbins).
!   nbins       -- integer: Beamlet size.
!   n_clamp     -- integer: Running count of nbl < 1 clamps, accumulated across calls.
!
! Output:
!   theta(:)    -- real(rp): Phases with the noise kicks added.
!   n_clamp     -- integer: Incremented once per clamped (harmonic, beamlet) draw.
!-

subroutine fel_fawley_noise (theta, weight, n, nbins, n_clamp)

real(rp) theta(:), weight(:)
integer n, nbins, n_clamp

real(rp), allocatable :: kick(:)
real(rp) nbl, u, phi, an
integer nharm, mbase, ih, ib, im, ip

!

nharm = (nbins - 1) / 2
mbase = n / nbins

allocate (kick(n))
kick = 0

do ih = 0, nharm - 1
  do ib = 1, mbase
    nbl = nbins * weight((ib-1)*nbins + 1) / e_charge
    if (nbl < 1) then
      nbl = 1
      n_clamp = n_clamp + 1
    endif
    call ran_uniform (u)
    phi = twopi * u
    call ran_uniform (u)
    an = sqrt(-log(u) / nbl) * 2 / real(ih+1, rp)
    if (an > twopi) an = mod(an, twopi)
    do im = 1, nbins
      ip = (ib-1)*nbins + im
      kick(ip) = kick(ip) - an * sin(theta(ip) * (ih+1) + phi)
    enddo
  enddo
enddo

theta(1:n) = theta(1:n) + kick(1:n)

end subroutine fel_fawley_noise

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_slice_to_bunch (beam, sl, ele, bunch, err_flag)
!
! Routine to convert a packed slice to a Bmad bunch_struct: plain copies, since the
! stored coordinates are coord_struct's. The element's p0c must match the beam's
! normalization. A mismatch is refused, not rescaled, since nothing on this path
! changes the reference momentum.
!
! Input:
!   beam        -- fel_beam_struct: The beam (for p0c and the chart assertion).
!   sl          -- fel_slice_struct: Slice to convert.
!   ele         -- ele_struct: Element at whose upstream end the coords are initialized.
!
! Output:
!   bunch       -- bunch_struct: The slice as a Bmad bunch.
!   err_flag    -- logical: Set True on error, False otherwise.
!-

subroutine fel_slice_to_bunch (beam, sl, ele, bunch, err_flag)

type (fel_beam_struct) beam
type (fel_slice_struct) sl
type (ele_struct) ele
type (bunch_struct) bunch
real(rp) vec(6)
integer ip
logical err_flag
character(*), parameter :: r_name = 'fel_slice_to_bunch'

!

err_flag = .true.

call fel_assert_averaged_chart (beam, 'fel_slice_to_bunch AT ELEMENT ' // trim(ele%name), err_flag)
if (err_flag) return
err_flag = .true.

if (abs(ele%value(p0c$) - beam%p0c) > 1e-10_rp * beam%p0c) then
  call out_io (s_error$, r_name, 'ELEMENT p0c \es20.12\ DOES NOT MATCH BEAM p0c \es20.12\ ', &
               'AT ELEMENT: ' // trim(ele%name), r_array = [ele%value(p0c$), beam%p0c])
  return
endif

if (allocated(bunch%particle)) then
  if (size(bunch%particle) /= sl%n) deallocate(bunch%particle)
endif
if (.not. allocated(bunch%particle)) allocate(bunch%particle(sl%n))

do ip = 1, sl%n
  vec = [sl%x(ip), sl%px(ip), sl%y(ip), sl%py(ip), sl%z(ip), sl%pz(ip)]

  ! init_coord derives beta, t and state consistently. shift_vec6 exists for elements
  ! whose reference momentum changes and must not touch vec(6) here.
  call init_coord (bunch%particle(ip), vec, ele, upstream_end$, electron$, shift_vec6 = .false.)
  bunch%particle(ip)%charge = sl%weight(ip)
enddo

bunch%n_live = sl%n
bunch%charge_live = sum(sl%weight(1:sl%n))

err_flag = .false.

end subroutine fel_slice_to_bunch

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_bunch_to_slice (bunch, ele, sl, err_flag)
!
! Routine to copy a tracked Bmad bunch back into the packed slice. Plain copies. The
! phase bookkeeping is one phi0 update per element, done by the caller
! (fel_phi0_advance), not here, because it is per beam and not per particle.
!
! Input:
!   bunch       -- bunch_struct: Tracked bunch. Every particle must be alive.
!   ele         -- ele_struct: Element just tracked through (for the error message).
!   sl          -- fel_slice_struct: Destination slice. The fill count n is unchanged.
!
! Output:
!   sl          -- fel_slice_struct: Slice with the tracked coordinates copied back.
!   err_flag    -- logical: Set True on error (a lost particle), False otherwise.
!-

subroutine fel_bunch_to_slice (bunch, ele, sl, err_flag)

type (bunch_struct) bunch
type (ele_struct) ele
type (fel_slice_struct) sl
integer ip
logical err_flag
character(*), parameter :: r_name = 'fel_bunch_to_slice'

!

err_flag = .true.

do ip = 1, sl%n
  if (bunch%particle(ip)%state /= alive$) then
    call out_io (s_error$, r_name, 'PARTICLE \i0\ LOST TRACKING THROUGH: ' // trim(ele%name), &
                                   i_array = [ip])
    return
  endif

  sl%x(ip)  = bunch%particle(ip)%vec(1)
  sl%px(ip) = bunch%particle(ip)%vec(2)
  sl%y(ip)  = bunch%particle(ip)%vec(3)
  sl%py(ip) = bunch%particle(ip)%vec(4)
  sl%z(ip)  = bunch%particle(ip)%vec(5)
  sl%pz(ip) = bunch%particle(ip)%vec(6)
  sl%weight(ip) = bunch%particle(ip)%charge
enddo

err_flag = .false.

end subroutine fel_bunch_to_slice

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_concat_slices (beam, ele, bunch, beta0, err_flag)
!
! Routine to concatenate all slices into ONE bunch_struct in GLOBAL window coordinates,
! so Bmad's whole-bunch machinery (the sr wake) sees the beam head to tail:
!
!   z_global = z_local + beta * (islice-1) * slice_spacing
!
! Higher slice index is the window HEAD (larger Bmad z). The direction and the formula
! are not new conventions. fel_migrate_slices documents "positive theta drift moves
! toward higher slice index, the head", and shifts a mover's z by exactly
! -atar*beta*slice_spacing (which makes z_global the migration INVARIANT).
! fel_wake_update's convolution collects current(is+i) from higher indices, the wake
! trailing its source. (An earlier specification guessed the opposite sign. The
! causality check and these two authorities pinned it.)
!
! beta0 records each particle's entry beta, slice-major in bunch order, for the exact
! inverse in fel_split_slices. Particles are stored slice-major (slice 1 first), and
! nothing downstream may permute bunch%particle storage. Bmad's wake routines order
! through the bunch%ix_z index array, never by moving particles, so the split can
! address segments by position.
!
! Input:
!   beam        -- fel_beam_struct: The beam to concatenate.
!   ele         -- ele_struct: Element at whose upstream end the coords are initialized.
!
! Output:
!   bunch       -- bunch_struct: All slices as one bunch, slice-major.
!   beta0(:)    -- real(rp), allocatable: Each particle's entry beta, in bunch order.
!   err_flag    -- logical: Set True on error, False otherwise.
!-

subroutine fel_concat_slices (beam, ele, bunch, beta0, err_flag)

type (fel_beam_struct), target :: beam
type (fel_slice_struct), pointer :: sl
type (ele_struct) ele
type (bunch_struct) bunch
real(rp), allocatable :: beta0(:)
real(rp) vec(6), beta, p0_mc
integer is, ip, k, ntot
logical err_flag
character(*), parameter :: r_name = 'fel_concat_slices'

!

err_flag = .true.

if (abs(ele%value(p0c$) - beam%p0c) > 1e-10_rp * beam%p0c) then
  call out_io (s_error$, r_name, 'ELEMENT p0c \es20.12\ DOES NOT MATCH BEAM p0c \es20.12\ ', &
               'AT ELEMENT: ' // trim(ele%name), r_array = [ele%value(p0c$), beam%p0c])
  return
endif

p0_mc = fel_p0_mc(beam)
ntot = 0
do is = 1, size(beam%slice)
  ntot = ntot + beam%slice(is)%n
enddo

if (allocated(bunch%particle)) then
  if (size(bunch%particle) /= ntot) deallocate(bunch%particle)
endif
if (.not. allocated(bunch%particle)) allocate(bunch%particle(ntot))
if (allocated(beta0)) then
  if (size(beta0) /= ntot) deallocate(beta0)
endif
if (.not. allocated(beta0)) allocate(beta0(ntot))

k = 0
do is = 1, size(beam%slice)
  sl => beam%slice(is)
  do ip = 1, sl%n
    k = k + 1
    beta = fel_beta_of(p0_mc, sl%pz(ip))
    beta0(k) = beta
    vec = [sl%x(ip), sl%px(ip), sl%y(ip), sl%py(ip), &
           sl%z(ip) + beta * (is-1) * beam%slice_spacing, sl%pz(ip)]
    call init_coord (bunch%particle(k), vec, ele, upstream_end$, electron$, shift_vec6 = .false.)
    bunch%particle(k)%charge = sl%weight(ip)
  enddo
enddo

bunch%n_live = ntot
bunch%charge_live = sum(bunch%particle(1:ntot)%charge)

err_flag = .false.

end subroutine fel_concat_slices

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_split_slices (bunch, ele, beam, beta0, hold_theta, err_flag)
!
! Routine to split a concatenated bunch back into the beam's slices: the exact inverse
! of fel_concat_slices. Each particle's stored entry offset beta0*(islice-1)*slice_spacing
! is subtracted from vec(5) and the segment copied back into its slice. Slice membership
! is by storage position (slice-major), which Bmad's wake routines never permute.
!
! hold_theta: rescale the restored z by beta_new/beta_old, so the derived phase
! theta = phi0 + ks*z/beta is untouched by an energy kick. This is Genesis's convention
! for wake energy loss (fel_wake_apply_slice does the same), used when the caller applied
! a pure kick with no transport (the in-wiggler wake). A track1_bunch passage uses
! hold_theta = false: there the element transported vec(5) in Bmad's own chart and the
! restored z IS the tracked z.
!
! Input:
!   bunch       -- bunch_struct: Whole-window bunch after the passage. All alive.
!   ele         -- ele_struct: Element just passed through (for the error message).
!   beam        -- fel_beam_struct: The beam whose slices receive the particles.
!   beta0(:)    -- real(rp), allocatable: Entry betas recorded by fel_concat_slices.
!   hold_theta  -- logical: Rescale z to hold theta fixed through an energy kick?
!
! Output:
!   beam        -- fel_beam_struct: Slices updated with the tracked coordinates.
!   err_flag    -- logical: Set True on error (a lost particle), False otherwise.
!-

subroutine fel_split_slices (bunch, ele, beam, beta0, hold_theta, err_flag)

type (bunch_struct) bunch
type (ele_struct) ele
type (fel_beam_struct), target :: beam
type (fel_slice_struct), pointer :: sl
real(rp), allocatable :: beta0(:)
real(rp) z_loc, beta_new, p0_mc
integer is, ip, k
logical hold_theta, err_flag
character(*), parameter :: r_name = 'fel_split_slices'

!

err_flag = .true.
p0_mc = fel_p0_mc(beam)

k = 0
do is = 1, size(beam%slice)
  sl => beam%slice(is)
  do ip = 1, sl%n
    k = k + 1
    if (bunch%particle(k)%state /= alive$) then
      call out_io (s_error$, r_name, 'PARTICLE \i0\ LOST THROUGH THE WHOLE-BUNCH PASSAGE OF: ' &
                   // trim(ele%name), i_array = [k])
      return
    endif

    z_loc = bunch%particle(k)%vec(5) - beta0(k) * (is-1) * beam%slice_spacing
    if (hold_theta) then
      beta_new = fel_beta_of(p0_mc, bunch%particle(k)%vec(6))
      z_loc = z_loc * beta_new / beta0(k)
    endif

    sl%x(ip)  = bunch%particle(k)%vec(1)
    sl%px(ip) = bunch%particle(k)%vec(2)
    sl%y(ip)  = bunch%particle(k)%vec(3)
    sl%py(ip) = bunch%particle(k)%vec(4)
    sl%z(ip)  = z_loc
    sl%pz(ip) = bunch%particle(k)%vec(6)
  enddo
enddo

err_flag = .false.

end subroutine fel_split_slices

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_migrate_slices (beam, ks, n_moved, charge_dropped, drop_re, drop_im, err_flag)
!
! Routine to move particles between slices when their ponderomotive phase leaves the
! slice window: the weighted generalization of Genesis's one4one-only Sorting::localSort
! (fel-physics.tex sec:migration). This port can offer it for ANY beam because each
! particle carries its own charge -- weights and migration are the same feature.
!
! The criterion is Genesis's, in this code's chart: the derived phase
! theta = phi0 + ks*z/beta is Genesis's bounded theta, the slice window is
! [0, slen) with slen = ks*slice_spacing = 2*pi*sample, and
! atar = floor(theta/slen) is the relative destination (localSort's exact formula:
! positive theta drift moves toward higher slice index, the head of the window). A
! mover's z shifts by exactly -atar*beta*slice_spacing, which changes theta by
! -atar*2*pi*sample. For integer sample the phase seen by every deposition and
! diagnostic is continuous across the move to rounding -- no wrap protocol, the
! unwrapped-z coordinate choice (module header) paying off. sample is asserted integer
! for that reason.
!
! Removal is Genesis's swap-with-last, and the while-loop re-examines the swapped-in
! particle exactly as localSort does. A particle appended to a HIGHER slice is
! re-examined when the scan reaches that slice (its adjusted theta then lies inside the
! window, so it stays). One appended to a LOWER slice waits for the next call, as in
! Genesis. Particles whose destination lies beyond the window are DROPPED WITH THEIR
! CHARGE COUNTED into charge_dropped. Genesis discards them silently at the world
! edges (sec:migration). The accounting is this port's deviation, chosen so
! conservation is checkable.
!
! Serial by design: called between the parallel regions at the caller's stride, so
! thread-count independence is untouched. Single-slice beams return immediately (a
! steady-state slice is periodic, so migration has no meaning).
!
! Input:
!   beam        -- fel_beam_struct: The beam.
!   ks          -- real(rp): Radiation wavenumber twopi/wavelength [1/m].
!
! Output:
!   beam             -- fel_beam_struct: Particles re-sliced. Slice fill counts updated.
!   n_moved          -- integer: Number of particles moved between slices (this call).
!   charge_dropped   -- real(rp): Charge of particles dropped off the window ends [C].
!   drop_re, drop_im -- real(rp): Weighted phasor sum(w*e^{i theta}) of the dropped
!                        particles, at their phase when dropped. Lets the caller verify
!                        exact phase continuity including drops:
!                        S_before = S_after + S_dropped to rounding.
!   err_flag         -- logical: Set True on error, False otherwise.
!-

subroutine fel_migrate_slices (beam, ks, n_moved, charge_dropped, drop_re, drop_im, err_flag)

type (fel_beam_struct), target :: beam
type (fel_slice_struct), pointer :: sl, sd
real(rp) ks, charge_dropped, drop_re, drop_im
integer n_moved
logical err_flag

real(rp) p0_mc, slen, sample, theta, beta, z_new
integer ia, ib, il, nslice, atar, idest
character(*), parameter :: r_name = 'fel_migrate_slices'

!

err_flag = .true.
n_moved = 0
charge_dropped = 0
drop_re = 0
drop_im = 0

nslice = size(beam%slice)
if (nslice < 2) then
  err_flag = .false.
  return
endif

sample = beam%slice_spacing / beam%wavelength
if (abs(sample - nint(sample)) > 1e-9_rp * sample) then
  call out_io (s_error$, r_name, &
        'MIGRATION NEEDS AN INTEGER sample (SLICE SPACING OVER WAVELENGTH): PHASE CONTINUITY', &
        'ACROSS A MOVE HOLDS ONLY THEN. GOT: \es16.8\ ', r_array = [sample])
  return
endif

p0_mc = fel_p0_mc(beam)
slen = ks * beam%slice_spacing            ! Window length in phase: 2*pi*sample.

do ia = 1, nslice
  sl => beam%slice(ia)
  ib = 1
  do while (ib <= sl%n)
    beta = fel_beta_of(p0_mc, sl%pz(ib))
    theta = beam%phi0 + ks * sl%z(ib) / beta
    atar = int(floor(theta / slen))

    if (atar == 0) then
      ib = ib + 1
      cycle
    endif

    idest = ia + atar

    if (idest < 1 .or. idest > nslice) then
      charge_dropped = charge_dropped + sl%weight(ib)
      drop_re = drop_re + sl%weight(ib) * cos(theta)
      drop_im = drop_im + sl%weight(ib) * sin(theta)
    else
      sd => beam%slice(idest)
      if (sd%n + 1 > size(sd%x)) call fel_slice_reallocate (sd, max(sd%n + 1, (3 * sd%n) / 2))
      sd%n = sd%n + 1
      z_new = sl%z(ib) - atar * beta * beam%slice_spacing
      sd%x(sd%n) = sl%x(ib);   sd%px(sd%n) = sl%px(ib)
      sd%y(sd%n) = sl%y(ib);   sd%py(sd%n) = sl%py(ib)
      sd%z(sd%n) = z_new;      sd%pz(sd%n) = sl%pz(ib)
      sd%weight(sd%n) = sl%weight(ib)
      n_moved = n_moved + 1
    endif

    ! Swap-with-last removal, re-examining the swapped-in particle (localSort's loop).

    il = sl%n
    sl%x(ib) = sl%x(il);   sl%px(ib) = sl%px(il)
    sl%y(ib) = sl%y(il);   sl%py(ib) = sl%py(il)
    sl%z(ib) = sl%z(il);   sl%pz(ib) = sl%pz(il)
    sl%weight(ib) = sl%weight(il)
    sl%n = sl%n - 1
  enddo
enddo

err_flag = .false.

end subroutine fel_migrate_slices

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Function fel_phi0_rate (ks, ku_like, p0_mc) result (rate)
!
! Routine to return dphi0/ds, the advance rate of the common ponderomotive reference
! phase:
!
!   dphi0/ds = ku_like + ks*(1 - 1/beta0)
!
! where ku_like is the undulator wavenumber inside an undulator, or Genesis's drift
! surrogate ks/(2*gamma0^2) in field-free regions, supplied by the caller
! (fel-physics.tex sec:chart), and beta0 is the reference beta from p0_mc. Combined with the
! per-particle theta_j = phi0 - ks*tau_j this reproduces Genesis's
! dtheta/ds = ks*(1 - 1/beta_z) + ku_like exactly, because
! -ks*dtau/ds = -ks*(1/beta_z - 1/beta0).
!
! The 1 - 1/beta0 factor is formed cancellation-free:
! 1 - 1/beta0 = -1/(beta0*(gamma0b^2)*(1+beta0)) with gamma0b^2 = p0_mc^2 + 1.
!
! Input:
!   ks       -- real(rp): Radiation wavenumber twopi/wavelength [1/m].
!   ku_like  -- real(rp): Undulator wavenumber, or the drift surrogate ks/(2*gamma0^2) [1/m].
!   p0_mc    -- real(rp): Reference momentum p0/(m_e c). See fel_p0_mc.
!
! Output:
!   rate     -- real(rp): dphi0/ds [rad/m].
!-

function fel_phi0_rate (ks, ku_like, p0_mc) result (rate)

real(rp) ks, ku_like, p0_mc, rate
real(rp) gamma0b_sq, beta0

!

gamma0b_sq = p0_mc**2 + 1
beta0 = p0_mc / sqrt(gamma0b_sq)
rate = ku_like - ks / (beta0 * gamma0b_sq * (1 + beta0))

end function fel_phi0_rate

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_slice_diag (beam, sl, ks, diag)
!
! Routine to compute the per-slice beam diagnostics, weighted. Genesis's DiagBeam::calc
! definitions with 1/N generalized to w_j/sum(w) (two-pass variances by REQUIREMENT:
! fel-physics.tex sec:numerics). Uniform weights reproduce Genesis exactly up to the
! regrouped arithmetic. Positions and sizes are of x, y directly. mean_px, mean_py are
! in Bmad's normalization P/p0 (multiply by fel_p0_mc for Genesis's gamma*beta). Adds
! N_eff and the derived current.
!
! Input:
!   beam    -- fel_beam_struct: The beam (for phi0, p0c, and slice_spacing).
!   sl      -- fel_slice_struct: Slice to analyze.
!   ks      -- real(rp): Radiation wavenumber twopi/wavelength [1/m].
!
! Output:
!   diag    -- fel_slice_diag_struct: The weighted diagnostics. Zeroed if the slice
!               carries no charge.
!-

subroutine fel_slice_diag (beam, sl, ks, diag)

type (fel_beam_struct) beam
type (fel_slice_struct) sl
type (fel_slice_diag_struct) diag
real(rp) ks
real(rp) g1, g2, x1, x2, y1, y2, px1, py1, br, bi, wsum, w2sum, w, gam, theta, beta, p_mc, p0_mc
integer ip

!

g1 = 0; g2 = 0; x1 = 0; x2 = 0; y1 = 0; y2 = 0; px1 = 0; py1 = 0
br = 0; bi = 0; wsum = 0; w2sum = 0
p0_mc = fel_p0_mc(beam)

do ip = 1, sl%n
  w = sl%weight(ip)
  p_mc = p0_mc * (1 + sl%pz(ip))         ! gamma and theta's beta = P/gamma share one
  gam = sqrt(p_mc**2 + 1)                ! sqrt. Bit-identical to fel_gamma_of/fel_theta.
  beta = p_mc / gam
  theta = beam%phi0 + ks * sl%z(ip) / beta

  wsum = wsum + w
  w2sum = w2sum + w*w
  x1 = x1 + w * sl%x(ip)
  y1 = y1 + w * sl%y(ip)
  g1 = g1 + w * gam
  px1 = px1 + w * sl%px(ip)
  py1 = py1 + w * sl%py(ip)
  br = br + w * cos(theta)
  bi = bi + w * sin(theta)
enddo

if (wsum <= 0) then
  diag = fel_slice_diag_struct()
  return
endif

g1 = g1/wsum
x1 = x1/wsum
y1 = y1/wsum

! Second pass for the variances (fel-physics.tex sec:numerics): the one-pass <v^2> - <v>^2 form loses
! sigma to cancellation once the mean is large against the spread. For gamma ~ 1e4
! with sigma ~ 1e-3 the one-pass noise is of order sigma itself, seen as spurious
! sigma_gamma jitter under uniform wake kicks before this was rewritten.

do ip = 1, sl%n
  w = sl%weight(ip)
  gam = fel_gamma_of(p0_mc, sl%pz(ip))
  g2 = g2 + w * (gam - g1)**2
  x2 = x2 + w * (sl%x(ip) - x1)**2
  y2 = y2 + w * (sl%y(ip) - y1)**2
enddo

diag%mean_energy = g1 * m_electron
diag%sigma_energy = sqrt(g2/wsum) * m_electron
diag%mean_x = x1
diag%sigma_x = sqrt(x2/wsum)
diag%mean_y = y1
diag%sigma_y = sqrt(y2/wsum)
diag%mean_px = px1/wsum
diag%mean_py = py1/wsum
diag%bunching = sqrt((br/wsum)**2 + (bi/wsum)**2)
diag%bunching_phase = atan2(bi/wsum, br/wsum)
diag%n_eff = wsum*wsum / w2sum
diag%current = c_light * wsum / beam%slice_spacing

end subroutine fel_slice_diag

end module fel_beam_mod
