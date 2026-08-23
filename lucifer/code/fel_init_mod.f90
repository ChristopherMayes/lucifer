!+
! Module fel_init_mod
!
! The starting state of a run: fel_init_beam (a Genesis dump, an imported
! distribution, or the generated quiet start) and fel_init_wavefront (a field dump,
! an openPMD wavefront, or the generated Gaussian seed, plus the harmonic entries).
! Library contract: errors return through err_flag; nothing here stops. The print
! lines are unchanged from when this code lived in the driver.
!-

module fel_init_mod

use fel_struct
use fel_io_mod
use beam_mod

implicit none

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_init_beam (run, err_flag)
!
! Routine to build the beam: read a Genesis particle dump (beam_file), import a
! distribution (dist_file or use_beam_init: the resample of fel_import_mod), or generate
! the quiet start from beam_init. Applies the beam-side check instruments
! (split_weights, swap_beam_xy) and sets run%nslice. One seed (global%ran_seed)
! governs generation, resampling and noise, exactly as before the split. Errors
! return through err_flag; nothing here stops.
!
! Input:
!   run       -- fel_run_struct: Run state after fel_setup_lattice (needs %gamma0, %lat and
!                  the parsed beam inputs).
!
! Output:
!   run       -- fel_run_struct: Beam built (%fbeam, %nslice).
!   err_flag  -- logical: Set True if there is an error. False otherwise.
!-

subroutine fel_init_beam (run, err_flag)

type (fel_run_struct), target :: run
logical err_flag

type (lat_struct), pointer :: lat
type (branch_struct), pointer :: branch
type (fel_beam_struct), pointer :: fbeam
type (beam_init_struct), pointer :: beam_init
type (fel_import_param_struct) imp
character(400) beam_file, dist_file, write_dist_file, write_opmd_file
character(400) out_root
character(400) field_file(9)
logical use_beam_init, shotnoise, gen_test_weights, imp_split_weights
logical split_weights, swap_beam_xy, err
integer nbins, ran_seed, is, ih
real(rp) gamma0, lambda0, window_length, seed_power, seed_waist_size, grid_half_width
integer window_sample, grid_n_pts, npart_gen
real(rp) delgam_gen, tw_beta_x, tw_alpha_x, tw_beta_y, tw_alpha_y
character(*), parameter :: r_name = 'fel_init_beam'

!

err_flag = .false.
lat => run%lat
branch => lat%branch(0)
fbeam => run%fbeam
beam_init => run%beam_init
imp = run%imp
gamma0 = run%gamma0
beam_file = run%bparam%beam_file
dist_file = run%bparam%dist_file
write_dist_file = run%bparam%write_dist_file
write_opmd_file = run%bparam%write_opmd_file
use_beam_init = run%bparam%use_beam_init
nbins = run%bparam%nbins
shotnoise = run%bparam%shotnoise
gen_test_weights = run%bparam%gen_test_weights
imp_split_weights = run%bparam%imp_split_weights
split_weights = run%bparam%split_weights
swap_beam_xy = run%bparam%swap_beam_xy
ran_seed = run%global%ran_seed
out_root = run%global%out_root
field_file = run%field_file
lambda0 = run%winit%lambda0
window_length = run%winit%window_length
window_sample = run%winit%window_sample
grid_n_pts = run%winit%grid_n_pts
grid_half_width = run%winit%grid_half_width
seed_power = run%winit%seed_power
seed_waist_size = run%winit%seed_waist_size

! ONE reference energy, and the lattice is it: gamma0 = e_tot/m_e c^2 from the lattice
! header, never a namelist input. There used to be a namelist gamma0 for Genesis-deck
! symmetry; the first external user fed it a hand-rounded value against a round lattice
! e_tot, the two disagreed at 1.4e-9, and the run died mid-tracking on the seam's
! backstop p0c check with raw numbers -- the FEL physics ran on one reference while
! Bmad's momenta were normalized by the other. Two specifications of one truth is the
! defect; the redundant one was removed (parameters live on the lattice).

if ((beam_file == '') .neqv. (field_file(1) == '')) then
  call out_io (s_error$, r_name, 'GIVE BOTH BEAM_FILE AND FIELD_FILE, OR NEITHER (TO GENERATE).')
  err_flag = .true.;  return
endif
if (field_file(1) == '' .and. any(field_file(2:) /= '')) then
  call out_io (s_error$, r_name, 'HARMONIC FIELD FILES NEED THE FUNDAMENTAL IN FIELD_FILE(1).')
  err_flag = .true.;  return
endif
if (beam_file /= '' .and. (dist_file /= '' .or. use_beam_init)) then
  call out_io (s_error$, r_name, 'DUMP FILES AND A DISTRIBUTION IMPORT ARE MUTUALLY EXCLUSIVE.')
  err_flag = .true.;  return
endif
if (dist_file /= '' .and. use_beam_init) then
  call out_io (s_error$, r_name, 'GIVE DIST_FILE OR USE_BEAM_INIT, NOT BOTH.')
  err_flag = .true.;  return
endif

if (beam_file /= '') then
  call fel_read_genesis4_beam (fbeam, beam_file, gamma0, err)
  if (err) then
    err_flag = .true.;  return
  endif
elseif (dist_file /= '' .or. use_beam_init) then
  call import_initial_state ()
  if (err_flag) return
else
  call generate_initial_state ()
  if (err_flag) return
endif

if (split_weights) call do_split_weights (fbeam)
run%nslice = size(fbeam%slice)

! Check instrument (check_two_polarization.py's rotation identity): swap the beam's
! transverse planes after generation, (x,px) <-> (y,py).

if (swap_beam_xy) then
  do is = 1, run%nslice
    call swap_arrays (fbeam%slice(is)%x, fbeam%slice(is)%y)
    call swap_arrays (fbeam%slice(is)%px, fbeam%slice(is)%py)
  enddo
endif

!------------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!+
! Subroutine generate_initial_state ()
!
! Routine to generate the quiet-start beam from the beam_init description (manual
! sec:loading): matched Gaussian transverse planes on the lattice Twiss, beamlets on a
! uniform ponderomotive phase grid, the derived per-slice current, and optional Fawley
! shot noise.
!-

subroutine generate_initial_state ()

type (fel_slice_struct), pointer :: sl
real(rp) p0_mc, ks_l, eg_x, eg_y, u, v, x, xp, y, yp, gam, p_mc, beta, pz, theta0
real(rp) dx_grid, w_part, e0, xg, yg, wsum, w2sum, n_lambda, n_eff, floor_b2, target_b2
real(rp) phi, an, nbl, br, bi
real(rp), allocatable :: theta_work(:), beta_work(:), kick(:), cur_gen(:)
real(rp) nl_min, nl_max, neff_min, neff_max, floor_max
real(rp) spacing_gen, zlen_gen, s_i
integer ib, im, ip, mbase, ix, iy, is_g, nslice_gen, ih, nharm, n_clamp
logical flat_z
character(*), parameter :: r_name = 'generate_initial_state'

!

if (lambda0 <= 0) then
  call out_io (s_error$, r_name, 'GENERATION NEEDS LAMBDA0 > 0 (REQUIRED BY DECISION; NOT DEFAULTED', &
                                 'FROM THE LATTICE RESONANCE -- THE FIRST UNDULATOR MAY BE OFF).')
  err_flag = .true.;  return
endif

call check_beam_init_contract ()

npart_gen = beam_init%n_particle
if (npart_gen < 1 .or. nbins < 1 .or. mod(npart_gen, nbins) /= 0) then
  call out_io (s_error$, r_name, 'BEAM_INIT%N_PARTICLE (MACROPARTICLES PER SLICE HERE) MUST BE A', &
                                 'POSITIVE MULTIPLE OF NBINS.')
  err_flag = .true.;  return
endif
if (beam_init%a_norm_emit <= 0 .or. beam_init%b_norm_emit <= 0) then
  call out_io (s_error$, r_name, 'BEAM_INIT%A_NORM_EMIT AND %B_NORM_EMIT MUST BE POSITIVE.')
  err_flag = .true.;  return
endif
if (beam_init%bunch_charge <= 0) then
  call out_io (s_error$, r_name, 'THE CURRENT IS DERIVED FROM THE DESCRIPTION; BEAM_INIT%BUNCH_CHARGE', &
                                 'MUST BE POSITIVE (THERE IS NO CURRENT PARAMETER).')
  err_flag = .true.;  return
endif

! The Twiss is the LATTICE's (one specification of one truth, as with e_tot and the
! import path's init_beam_distribution): read the beginning element, refuse by name
! when a lattice carries none.

tw_beta_x = branch%ele(0)%a%beta;  tw_alpha_x = branch%ele(0)%a%alpha
tw_beta_y = branch%ele(0)%b%beta;  tw_alpha_y = branch%ele(0)%b%alpha
if (tw_beta_x <= 0 .or. tw_beta_y <= 0) then
  call out_io (s_error$, r_name, 'THE LATTICE CARRIES NO BEGINNING TWISS (BEGINNING[BETA_A], ETC.);', &
                                 'THE GENERATED QUIET START IS MATCHED TO THE LATTICE, SO THE LATTICE MUST SAY.')
  err_flag = .true.;  return
endif
if (beam_init%sig_pz < 0 .or. seed_power < 0 .or. grid_n_pts < 3 .or. grid_half_width <= 0 .or. &
    window_sample < 1) then
  call out_io (s_error$, r_name, 'CHECK BEAM_INIT%SIG_PZ, SEED_POWER, GRID_N_PTS, GRID_HALF_WIDTH,', &
                                 'WINDOW_SAMPLE.')
  err_flag = .true.;  return
endif
if (seed_power > 0 .and. seed_waist_size <= 0) then
  call out_io (s_error$, r_name, 'SEED_WAIST_SIZE MUST BE POSITIVE WHEN SEED_POWER > 0.')
  err_flag = .true.;  return
endif

! The window and the per-slice current derive from the beam_init description (manual
! sec:loading): one bulk bunch, evaluated analytically at the slice centers. The
! default window covers the described bunch (as the import derives its window from
! real particles); window_length overrides it for slippage headroom and warns when it
! clips the bunch. sig_z = 0 is the steady state -- the whole charge in one slice
! window -- and is refused by name for time-dependent windows.

flat_z = .false.
select case (trim(beam_init%distribution_type(3)))
case ('', 'RAN_GAUSS', 'ran_gauss', 'Ran_Gauss')
case ('GRID', 'grid', 'Grid')
  flat_z = .true.
case default
  call out_io (s_error$, r_name, 'BEAM_INIT%DISTRIBUTION_TYPE(3) MUST BE RAN_GAUSS (GAUSSIAN BUNCH)', &
               'OR GRID (FLAT, BMAD''S UNIFORM) FOR THE QUIET-START GENERATOR, GOT: ' // &
               trim(beam_init%distribution_type(3)))
  err_flag = .true.;  return
end select

spacing_gen = window_sample * lambda0

if (flat_z) then
  zlen_gen = beam_init%grid(3)%x_max - beam_init%grid(3)%x_min
  if (zlen_gen <= 0) then
    call out_io (s_error$, r_name, 'A GRID (FLAT) Z-PLANE NEEDS BEAM_INIT%GRID(3)%X_MIN < %X_MAX.')
    err_flag = .true.;  return
  endif
elseif (beam_init%sig_z > 0) then
  zlen_gen = 8 * beam_init%sig_z          ! +-4 sigma covers the described bunch.
else
  zlen_gen = 0                            ! Steady state.
endif

if (window_length > 0) then
  if (zlen_gen == 0 .and. window_length > 1.5_rp * spacing_gen) then
    call out_io (s_error$, r_name, 'SIG_Z = 0 (THE STEADY-STATE DESCRIPTION) IS INVALID FOR A', &
                 'TIME-DEPENDENT WINDOW. GIVE THE BUNCH A LENGTH (SIG_Z, OR A GRID EXTENT).')
    err_flag = .true.;  return
  endif
  if (window_length < zlen_gen) then
    call out_io (s_warn$, r_name, 'window_length = \es10.3\ m CLIPS the described bunch (\es10.3\ m).', &
                 r_array = [window_length, zlen_gen])
  endif
  nslice_gen = max(1, nint(window_length / spacing_gen))
else
  nslice_gen = max(1, nint(zlen_gen / spacing_gen))
endif

if (shotnoise .and. nslice_gen < 2) then
  call out_io (s_error$, r_name, 'SHOTNOISE NEEDS A TIME-DEPENDENT WINDOW, THE SAME RULE AS GENESIS.')
  err_flag = .true.;  return
endif

mbase = npart_gen / nbins
if (gen_test_weights .and. mod(mbase, 2) /= 0) then
  call out_io (s_error$, r_name, 'GEN_TEST_WEIGHTS NEEDS AN EVEN NUMBER OF BEAMLETS.')
  err_flag = .true.;  return
endif

p0_mc = sqrt(gamma0**2 - 1)
delgam_gen = (p0_mc**2 / gamma0) * beam_init%sig_pz    ! beta0*p0_mc*sig_pz.

fbeam%p0c = p0_mc * m_electron
fbeam%phi0 = 0
fbeam%wavelength = lambda0
fbeam%slice_spacing = spacing_gen
fbeam%s0 = 0
fbeam%nbins = nbins
fbeam%one4one = .false.

if (allocated(fbeam%slice)) deallocate(fbeam%slice)
allocate (fbeam%slice(nslice_gen), cur_gen(nslice_gen))

! The derived per-slice current: flat Q*c/extent inside the GRID extent; Gaussian
! profile at the slice centers, bunch centered in the window; steady state = the
! whole charge in the one slice window, I = Q*c/spacing.

if (flat_z) then
  cur_gen = 0
  do is_g = 1, nslice_gen
    s_i = (is_g - 1) * spacing_gen - (nslice_gen - 1) * spacing_gen / 2
    if (abs(s_i) <= zlen_gen / 2) cur_gen(is_g) = beam_init%bunch_charge * c_light / zlen_gen
  enddo
elseif (zlen_gen > 0) then
  do is_g = 1, nslice_gen
    s_i = (is_g - 1) * spacing_gen - (nslice_gen - 1) * spacing_gen / 2
    cur_gen(is_g) = beam_init%bunch_charge * c_light / (sqrt(twopi) * beam_init%sig_z) * &
                    exp(-s_i**2 / (2 * beam_init%sig_z**2))
  enddo
else
  cur_gen(1) = beam_init%bunch_charge * c_light / spacing_gen
endif

call ran_seed_put (ran_seed)

ks_l = twopi / lambda0
eg_x = beam_init%a_norm_emit / p0_mc  ! Normalized emittance to geometric.
eg_y = beam_init%b_norm_emit / p0_mc

allocate (theta_work(npart_gen), beta_work(npart_gen))
n_clamp = 0
nl_min = huge(1.0_rp); nl_max = 0; neff_min = huge(1.0_rp); neff_max = 0; floor_max = 0

do is_g = 1, nslice_gen
  sl => fbeam%slice(is_g)
  call fel_slice_reallocate (sl, npart_gen)
  sl%n = npart_gen
  w_part = cur_gen(is_g) * fbeam%slice_spacing / (c_light * npart_gen)

  ! Quiet start: mbase base samples, each replicated at nbins equally spaced
  ! ponderomotive phases (theta0 spread on a uniform grid within one beamlet spacing),
  ! so bunching harmonics below nbins vanish to roundoff. Weights and coordinates
  ! follow fel_read_genesis4_beam: z = beta*theta/ks with phi0 = 0,
  ! weight = I*slice_spacing/(c*npart). theta and beta are held in work arrays so noise
  ! can kick the phases before the z conversion.

  ip = 0
  do ib = 1, mbase
    call ran_gauss (u);  call ran_gauss (v)
    x  = sqrt(eg_x * tw_beta_x) * u
    xp = sqrt(eg_x / tw_beta_x) * (v - tw_alpha_x * u)
    call ran_gauss (u);  call ran_gauss (v)
    y  = sqrt(eg_y * tw_beta_y) * u
    yp = sqrt(eg_y / tw_beta_y) * (v - tw_alpha_y * u)

    call ran_gauss (u)
    gam = gamma0 + delgam_gen * u
    p_mc = sqrt(gam**2 - 1)
    beta = p_mc / gam
    pz = (p_mc - p0_mc) / p0_mc

    theta0 = (ib - 0.5_rp) * twopi / (nbins * mbase)

    do im = 0, nbins - 1
      ip = ip + 1
      theta_work(ip) = theta0 + im * twopi / nbins
      beta_work(ip) = beta
      sl%x(ip) = x;   sl%px(ip) = xp
      sl%y(ip) = y;   sl%py(ip) = yp
      sl%pz(ip) = pz
      sl%weight(ip) = w_part
    enddo
  enddo

  ! Validation knob: alternate beamlet weights 0.25x/1.75x, charge preserving, uniform
  ! within each beamlet so the quiet cancellation is untouched. Exercises every
  ! weighted-noise path (the asymmetry is strong enough that using a slice-uniform
  ! electron count where the beamlet's charge belongs mis-sets <|b|^2> by 56 percent,
  ! far outside the statistical check); not a physics input.

  if (gen_test_weights) then
    do ib = 1, mbase
      sl%weight((ib-1)*nbins+1 : ib*nbins) = &
              sl%weight((ib-1)*nbins+1 : ib*nbins) * (1 + 0.75_rp * (-1)**ib)
    enddo
  endif

  ! Noise bookkeeping (manual sec:noise): real electrons N_lambda = charge/e, effective
  ! macroparticle number N_eff = (sum w)^2/sum w^2, both per slice.

  wsum = sum(sl%weight(1:npart_gen))
  w2sum = sum(sl%weight(1:npart_gen)**2)
  n_lambda = wsum / e_charge
  n_eff = wsum**2 / w2sum
  nl_min = min(nl_min, n_lambda);  nl_max = max(nl_max, n_lambda)
  neff_min = min(neff_min, n_eff); neff_max = max(neff_max, n_eff)

  ! Zero-current slices (Gaussian tails, outside a flat extent) carry no noise --
  ! Genesis's own zero-current skip, shared with the import.

  if (shotnoise .and. wsum > 0) then

    ! The N_eff guard: measure the pre-noise quiet floor. A representation whose floor
    ! is not far below the target 1/N_lambda cannot carry physical noise -- imposing on
    ! top would give a silently wrong startup level. The sweep covers EVERY harmonic the
    ! beamlet structure can resolve (1..nbins-1), not just the imposed ones: an
    ! unquiet weight pattern can park its floor on a harmonic the imposition never
    ! touches (an alternating within-beamlet pattern lands exactly on nbins/2, found
    ! by the guard's own mutation test) and still corrupt the dynamics through the
    ! nonlinear phase evolution.

    target_b2 = 1 / n_lambda
    floor_b2 = 0
    do ih = 1, nbins - 1
      br = 0; bi = 0
      do ip = 1, npart_gen
        br = br + sl%weight(ip) * cos(ih * theta_work(ip))
        bi = bi + sl%weight(ip) * sin(ih * theta_work(ip))
      enddo
      floor_b2 = max(floor_b2, (br**2 + bi**2) / wsum**2)
    enddo
    floor_max = max(floor_max, floor_b2 * n_lambda)

    if (floor_b2 > 0.01_rp * target_b2) then
      call out_io (s_error$, r_name, 'SLICE \i0\ : THE QUIET-START FLOOR IS NOT FAR BELOW THE', &
                   'PHYSICAL SHOT-NOISE LEVEL -- THIS REPRESENTATION CANNOT CARRY THE REQUESTED NOISE.', &
                   'MAX_H |B(H)|^2 = \es10.2\ VS TARGET 1/N_LAMBDA = \es10.2\ ', &
                   'N_EFF = \es10.2\ N_LAMBDA = \es10.2\ ', &
                   i_array = [is_g], r_array = [floor_b2, target_b2, n_eff, n_lambda])
      err_flag = .true.;  return
    endif

    ! Fawley-style shot noise: fel_fawley_noise (fel_beam_mod), the ShotNoise
    ! transcription generalized to weights, shared with the distribution import so the
    ! two paths stay one implementation. Draw order is unchanged from when this block
    ! lived inline here (two ran_uniform per harmonic per beamlet, Genesis's loops).

    call fel_fawley_noise (theta_work(1:npart_gen), sl%weight(1:npart_gen), npart_gen, nbins, n_clamp)
  endif

  ! To the stored chart: z = beta*theta/ks with phi0 = 0, beta of the base sample.

  do ip = 1, npart_gen
    sl%z(ip) = beta_work(ip) * theta_work(ip) / ks_l
  enddo
enddo

deallocate (theta_work, beta_work)

if (shotnoise) then
  call out_io (s_info$, r_name, 'Shot noise imposed on \i0\ slices.', &
               '  N_lambda per slice: \es10.3\ to \es10.3\ ', &
               '  N_eff per slice:    \es10.3\ to \es10.3\ ', &
               '  worst quiet floor, |b|^2 * N_lambda: \es10.2\ ', &
               i_array = [nslice_gen], r_array = [nl_min, nl_max, neff_min, neff_max, floor_max])
  if (n_clamp > 0) then
    call out_io (s_warn$, r_name, '\i0\ beamlet draws had fewer than one real electron', &
                 '(nbl clamped to 1, as Genesis does silently). The noise level in those', &
                 'beamlets is not physical; use fewer macroparticles or more charge.', &
                 i_array = [n_clamp])
  endif
endif

end subroutine generate_initial_state

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine check_beam_init_contract ()
!
! Routine to check the beam_init contract: the quiet-start generator honors the beam_init
! fields in the header table and REFUSES BY NAME every other field that is set -- a
! standard structure that silently dropped fields would be worse than a custom one. (The
! import path is exempt: init_beam_distribution honors everything Bmad honors.)
! renorm_center/renorm_sigma, random_engine defaults and n_bunch = 0/1 are generation
! details with no analytic counterpart and are accepted at their defaults only.
!-

subroutine check_beam_init_contract ()

character(60) bad
character(*), parameter :: r_name = 'check_beam_init_contract'

!

bad = ''
if (beam_init%position_file /= '')                          bad = 'position_file'
if (beam_init%a_emit /= 0 .or. beam_init%b_emit /= 0)       bad = 'a_emit/b_emit (use a_norm_emit/b_norm_emit)'
if (beam_init%dPz_dz /= 0)                                  bad = 'dPz_dz'
if (any(beam_init%center /= 0))                             bad = 'center'
if (any(beam_init%spin /= 0))                               bad = 'spin'
if (any(beam_init%center_jitter /= 0))                      bad = 'center_jitter'
if (any(beam_init%emit_jitter /= 0))                        bad = 'emit_jitter'
if (beam_init%sig_z_jitter /= 0)                            bad = 'sig_z_jitter'
if (beam_init%sig_pz_jitter /= 0)                           bad = 'sig_pz_jitter'
if (beam_init%t_offset /= 0)                                bad = 't_offset'
if (beam_init%dt_bunch /= 0)                                bad = 'dt_bunch'
if (beam_init%n_bunch > 1)                                  bad = 'n_bunch'
if (beam_init%ix_turn /= 0)                                 bad = 'ix_turn'
if (beam_init%full_6D_coupling_calc)                        bad = 'full_6D_coupling_calc'
if (beam_init%use_particle_start)                           bad = 'use_particle_start'
if (beam_init%use_t_coords)                                 bad = 'use_t_coords'
if (beam_init%file_name /= '')                              bad = 'file_name'
if (beam_init%random_engine /= '' .and. beam_init%random_engine /= 'pseudo') bad = 'random_engine'
if (beam_init%random_gauss_converter /= '' .and. beam_init%random_gauss_converter /= 'ziggurat') &
                                                            bad = 'random_gauss_converter'
if (beam_init%random_sigma_cutoff /= -1)                    bad = 'random_sigma_cutoff'
if (beam_init%species /= '' .and. beam_init%species /= 'electron') bad = 'species (electron only)'
if (trim(beam_init%distribution_type(1)) /= '' .and. trim(beam_init%distribution_type(1)) /= 'RAN_GAUSS' &
    .and. trim(beam_init%distribution_type(1)) /= 'ran_gauss') bad = 'distribution_type(1) (transverse: RAN_GAUSS only)'
if (trim(beam_init%distribution_type(2)) /= '' .and. trim(beam_init%distribution_type(2)) /= 'RAN_GAUSS' &
    .and. trim(beam_init%distribution_type(2)) /= 'ran_gauss') bad = 'distribution_type(2) (transverse: RAN_GAUSS only)'

if (bad /= '') then
  call out_io (s_error$, r_name, 'BEAM_INIT%' // trim(bad) // ' IS SET BUT NOT HONORED BY THE', &
                                 'QUIET-START GENERATOR (SEE THE HONORED-FIELDS TABLE IN THE PROGRAM HEADER).', &
                                 'REFUSING RATHER THAN SILENTLY IGNORING IT.')
  err_flag = .true.;  return
endif

end subroutine check_beam_init_contract

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine import_initial_state ()
!
! Routine to import a distribution (manual sec:import): a bunch_struct -- generated from
! Bmad's beam_init_struct (the native equivalent of Genesis's &beam description) or read
! from an openPMD-beamphysics file -- is resampled into FEL slices by the transcribed
! Genesis importdistribution method (fel_import_mod, where the algorithm and its
! provenance are documented). The seed field comes from the same generator as the
! built-in loader. The RNG-free outputs the exactness checks read -- the analysis
! moments and the per-slice current profile -- are printed at full precision.
!-

subroutine import_initial_state ()

type (beam_struct), target :: beam_b
type (bunch_struct), pointer :: bp
real(rp) moments(11)
integer is_g, ip_g, n0, iu_i
logical err_i
character(400) line
character(*), parameter :: r_name = 'import_initial_state'

!

if (lambda0 <= 0) then
  call out_io (s_error$, r_name, 'IMPORT NEEDS LAMBDA0 > 0.')
  err_flag = .true.;  return
endif
if (window_sample < 1) then
  call out_io (s_error$, r_name, 'WINDOW_SAMPLE MUST BE A POSITIVE INTEGER (GENESIS''S SAMPLE).')
  err_flag = .true.;  return
endif

! One seed governs the whole import: the bunch generation, the resampler's draws and
! the shot noise. Seeding AFTER generation was the first mutation this path caught in
! development -- every run then imports a different bunch, and the split-weight and
! thread-determinism checks both fail on what looks like resampler noise.

call ran_seed_put (ran_seed)

if (use_beam_init) then
  if (beam_init%n_particle < 1) then
    call out_io (s_error$, r_name, 'BEAM_INIT%N_PARTICLE MUST BE POSITIVE.')
    err_flag = .true.;  return
  endif
  beam_init%n_bunch = 1
  call init_beam_distribution (branch%ele(0), lat%param, beam_init, beam_b, err_i)
  if (err_i) then
    err_flag = .true.;  return
  endif
  call out_io (s_info$, r_name, 'Generated \i0\ particles from beam_init.', &
               i_array = [size(beam_b%bunch(1)%particle)])
else
  call hdf5_read_beam (dist_file, beam_b, err_i, branch%ele(0))
  if (err_i) then
    err_flag = .true.;  return
  endif
  call out_io (s_info$, r_name, 'Read \i0\ particles from: ' // trim(dist_file), &
               i_array = [size(beam_b%bunch(1)%particle)])
endif

bp => beam_b%bunch(1)

! Check knob: coincident split-weight copies before anything downstream sees the bunch.
! The current profile (weighted sums) and the analysis moments (unweighted, over
! coincident copies) must then be bit-identical to the unsplit run.

if (imp_split_weights) then
  n0 = size(bp%particle)
  call reallocate_bunch (bp, 2*n0, save = .true.)
  do ip_g = 1, n0
    bp%particle(n0+ip_g) = bp%particle(ip_g)
    bp%particle(n0+ip_g)%charge = 2 * bp%particle(ip_g)%charge / 3
    bp%particle(ip_g)%charge = bp%particle(ip_g)%charge / 3
  enddo
endif

if (write_dist_file /= '') then
  call fel_write_genesis4_distribution (bp, write_dist_file, err_i)
  if (err_i) then
    err_flag = .true.;  return
  endif
  call out_io (s_info$, r_name, 'Wrote Genesis distribution file: ' // trim(write_dist_file))
endif

if (write_opmd_file /= '') then
  call hdf5_write_beam (write_opmd_file, beam_b%bunch(1:1), .false., err_i, lat)
  if (err_i) then
    err_flag = .true.;  return
  endif
  call out_io (s_info$, r_name, 'Wrote openPMD-beamphysics file: ' // trim(write_opmd_file))
endif

! imp%npart and imp%nbins come from the imp block directly (the resample's own knobs;
! beam_init%n_particle is the BUNCH particle count on this path).
call fel_import_bunch (bp, gamma0, lambda0, window_sample * lambda0, imp, fbeam, err_i, moments)
if (err_i) then
  err_flag = .true.;  return
endif
call out_io (s_info$, r_name, 'Imported into \i0\ slices of \i0\ particles.', &
             i_array = [size(fbeam%slice), imp%npart])

! The RNG-free instruments the exactness checks read -- the analysis moments and the
! per-slice current profile -- go to a FILE, not stdout: stdout is for humans and
! nslice current lines are not (manual sec:program). Full precision, one row per slice.
! Written here, at import time, because load_only stops before tracking.

open (newunit = iu_i, file = trim(out_root) // '.import.txt', action = 'write')
write (iu_i, '(a)') '# The distribution import, at full precision. Machine-readable; stdout is not.'
write (iu_i, '(a, i0, a, i0)') '# nslice = ', size(fbeam%slice), '   npart_per_slice = ', imp%npart
write (iu_i, '(a)') '# moments: gavg xavg pxavg yavg pyavg ex ey bx by ax ay'
write (iu_i, '(a, 11es24.15e3)') 'moments', moments
write (iu_i, '(a)') '#  slice            current [A]'
do is_g = 1, size(fbeam%slice)
  write (iu_i, '(a, i0, a, es24.15e3)') 'current ', is_g, ' ', &
        c_light * sum(fbeam%slice(is_g)%weight(1:fbeam%slice(is_g)%n)) / fbeam%slice_spacing
enddo
close (iu_i)
call out_io (s_info$, r_name, 'Wrote ' // trim(out_root) // '.import.txt (moments and the current profile).')

end subroutine import_initial_state

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine do_split_weights (beam)
!
! Routine to replace each particle by two coincident copies with weights w/3 and 2w/3.
! The order -- all first copies, then all second copies -- keeps the original particles'
! storage order, which keeps the RK4 arithmetic per copy identical to the unsplit run.
!-

subroutine do_split_weights (beam)

type (fel_beam_struct), target :: beam
type (fel_slice_struct), pointer :: sp
integer is, ip, n0

do is = 1, size(beam%slice)
  sp => beam%slice(is)
  n0 = sp%n
  call fel_slice_reallocate (sp, 2*n0)
  do ip = 1, n0
    sp%x(n0+ip) = sp%x(ip);  sp%px(n0+ip) = sp%px(ip)
    sp%y(n0+ip) = sp%y(ip);  sp%py(n0+ip) = sp%py(ip)
    sp%z(n0+ip) = sp%z(ip);  sp%pz(n0+ip) = sp%pz(ip)
    sp%weight(n0+ip) = 2 * sp%weight(ip) / 3
    sp%weight(ip) = sp%weight(ip) / 3
  enddo
  sp%n = 2*n0
enddo

end subroutine do_split_weights

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine swap_arrays (a, b)
!
! Routine to swap the contents of the two arrays a and b.
!-

subroutine swap_arrays (a, b)

real(rp) a(:), b(:), tmp(size(a))

tmp = a;  a = b;  b = tmp

end subroutine swap_arrays

end subroutine fel_init_beam

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_init_wavefront (run, err_flag)
!
! Routine to build the field set: read the fundamental (field_file(1): a Genesis dump or
! an openPMD EXT_Wavefront, auto-detected by signature) or generate the Gaussian seed
! (wavefront_init; seed_power = 0 is a dark start); initialize the harmonic entries
! on the fundamental's grid; fill any from openPMD imports matched by photon energy;
! and check the beam/field window consistency. Needs the beam (fel_init_beam first).
! Errors return through err_flag; nothing here stops.
!
! Input:
!   run       -- fel_run_struct: Run state after fel_setup_lattice and fel_init_beam (needs
!                  %fbeam, %nslice, %two_pol, %n_harm and the parsed field inputs).
!
! Output:
!   run       -- fel_run_struct: Field set filled (%ffield(:)%wf, %ks).
!   err_flag  -- logical: Set True if there is an error. False otherwise.
!-

subroutine fel_init_wavefront (run, err_flag)

type (fel_run_struct), target :: run
logical err_flag

type (fel_field_struct), pointer :: ffield(:)
type (fel_beam_struct), pointer :: fbeam
type (wavefront_struct), pointer :: wf
character(400) field_file(9)
character(1) seed_polarization
logical two_pol, err
integer n_harm, nslice, ih, is
real(rp) lambda0, seed_power, seed_waist_size, grid_half_width
integer grid_n_pts
character(*), parameter :: r_name = 'fel_init_wavefront'

!

err_flag = .false.
ffield => run%ffield
fbeam => run%fbeam
wf => run%ffield(1)%wf
field_file = run%field_file
two_pol = run%two_pol
n_harm = run%n_harm
nslice = run%nslice
seed_polarization = run%winit%seed_polarization
lambda0 = run%winit%lambda0
seed_power = run%winit%seed_power
seed_waist_size = run%winit%seed_waist_size
grid_n_pts = run%winit%grid_n_pts
grid_half_width = run%winit%grid_half_width

if (field_file(1) /= '') then
  if (wavefront_file_is_openpmd(field_file(1))) then
    call fel_read_openpmd_into_field (run, field_file(1), 1, err)
    if (err) then
      err_flag = .true.;  return
    endif
  else
    call wavefront_read_genesis4 (wf, field_file(1), err)
    if (err) then
      err_flag = .true.;  return
    endif
  endif
  if (two_pol .and. .not. allocated(wf%Ey)) then
    allocate (wf%Ey(size(wf%Ex,1), size(wf%Ex,2), size(wf%Ex,3)))
    wf%Ey = 0
  endif
else
  call generate_seed_field (nslice)
  if (err_flag) return
endif

! Harmonic fields: the fundamental's grid and window, its wavelength / h, dark. An
! openPMD import (field_file entries 2+) fills the entry whose photon energy it
! carries; the fundamental's grid must match, per the same one-window rule the
! fundamental import obeys.

do ih = 2, n_harm
  call wavefront_init (ffield(ih)%wf, size(wf%Ex,1), size(wf%Ex,2), size(wf%Ex,3), &
                       wf%dx, wf%dy, wf%dz, wf%wavelength / ffield(ih)%harm, 'x', wf%ref_position)
enddo

do is = 2, 9
  if (field_file(is) == '') cycle
  if (.not. wavefront_file_is_openpmd(field_file(is))) then
    call out_io (s_error$, r_name, 'HARMONIC FIELD FILES MUST BE openPMD EXT_WAVEFRONT', &
                 '(THE GENESIS FORMAT CARRIES NO PHOTON ENERGY TO MATCH ON): ' // trim(field_file(is)))
    err_flag = .true.;  return
  endif
  call fel_import_harmonic_field (run, field_file(is), err)
  if (err) then
    err_flag = .true.;  return
  endif
enddo

run%ks = twopi / wf%wavelength

! The beam and field must describe the same time window: one field slice per beam slice,
! at the same wavelength. Checked, never assumed.

if (size(wf%Ex, 3) /= nslice) then
  call out_io (s_error$, r_name, 'BEAM HAS \i0\ SLICES BUT THE FIELD HAS \i0\ ', &
               i_array = [nslice, size(wf%Ex, 3)])
  err_flag = .true.;  return
endif
if (abs(wf%wavelength - fbeam%wavelength) > 1e-12_rp * fbeam%wavelength) then
  call out_io (s_error$, r_name, 'BEAM AND FIELD DISAGREE ON THE WAVELENGTH: \2es20.12\ ', &
               r_array = [fbeam%wavelength, wf%wavelength])
  err_flag = .true.;  return
endif

!------------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!+
! Subroutine generate_seed_field (nslice_f)
!
! Routine to generate the seed field: a Gaussian seed at its waist in every slice,
! E = E0*exp(-r^2/w0^2), intensity 1/e^2 radius w0, integrating to seed_power;
! seed_power = 0 is a dark start. Grid convention matches Genesis's dgrid: ngrid points
! spanning +-dgrid, dx = 2*dgrid/(ngrid-1), center on axis. Shared by the built-in
! generator and the distribution import (both make their own beam, neither brings a field).
!-

subroutine generate_seed_field (nslice_f)

integer nslice_f, ix, iy, is_g
real(rp) dx_grid, e0, xg, yg

!

if (grid_n_pts < 3 .or. grid_half_width <= 0) then
  call out_io (s_error$, r_name, 'CHECK GRID_N_PTS AND GRID_HALF_WIDTH.')
  err_flag = .true.;  return
endif
if (seed_power > 0 .and. seed_waist_size <= 0) then
  call out_io (s_error$, r_name, 'SEED_WAIST_SIZE MUST BE POSITIVE WHEN SEED_POWER > 0.')
  err_flag = .true.;  return
endif

dx_grid = 2 * grid_half_width / (grid_n_pts - 1)
call wavefront_init (wf, grid_n_pts, grid_n_pts, nslice_f, dx_grid, dx_grid, &
                     fbeam%slice_spacing, lambda0, 'x', 0.0_rp)

if (two_pol .and. .not. allocated(wf%Ey)) then
  allocate (wf%Ey(grid_n_pts, grid_n_pts, nslice_f))
  wf%Ey = 0
endif

if (seed_power > 0 .and. seed_polarization == 'y') then
  e0 = sqrt(4 * (mu_0_vac * c_light) * seed_power / (pi * seed_waist_size**2))
  do iy = 1, grid_n_pts
    yg = (iy - 1) * dx_grid - grid_half_width
    do ix = 1, grid_n_pts
      xg = (ix - 1) * dx_grid - grid_half_width
      wf%Ey(ix, iy, 1) = e0 * exp(-(xg**2 + yg**2) / seed_waist_size**2)
    enddo
  enddo
  do is_g = 2, nslice_f
    wf%Ey(:, :, is_g) = wf%Ey(:, :, 1)
  enddo
elseif (seed_power > 0) then
  e0 = sqrt(4 * (mu_0_vac * c_light) * seed_power / (pi * seed_waist_size**2))
  do iy = 1, grid_n_pts
    yg = (iy - 1) * dx_grid - grid_half_width
    do ix = 1, grid_n_pts
      xg = (ix - 1) * dx_grid - grid_half_width
      wf%Ex(ix, iy, 1) = e0 * exp(-(xg**2 + yg**2) / seed_waist_size**2)
    enddo
  enddo
  do is_g = 2, nslice_f
    wf%Ex(:, :, is_g) = wf%Ex(:, :, 1)
  enddo
endif

end subroutine generate_seed_field

end subroutine fel_init_wavefront

end module fel_init_mod
