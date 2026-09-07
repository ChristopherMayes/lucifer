!+
! Module fel_init_mod
!
! The starting state of a run: fel_init_beam (an openPMD dump of slices, or a beam_init
! bunch loaded by load_mode, sampled or kept) and fel_init_wavefront (an openPMD wavefront
! or the generated Gaussian seed, plus the harmonic entries).
! Library contract: errors return through err_flag, and nothing here stops. The print
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
! Routine to build the beam: read an openPMD particle dump of slices (beam_file), or load
! a beam_init bunch by load_mode (fel-physics.md sec-loading). In the sample mode a
! described bunch takes the analytic generator and a file, beam_init%position_file,
! takes the resampler of fel_import_mod. In the keep mode every particle of the bunch
! stays as a beamlet of copies. Applies the beam-side check instruments (split_weights,
! swap_beam_xy) and sets run%nslice. One seed (global%ran_seed) governs generation,
! resampling and noise. Errors return through err_flag, and nothing here stops.
!
! Input:
!   run       -- fel_run_struct: Run state after fel_setup_lattice (needs %gamma0, %lat and
!                  the parsed beam inputs).
!
! Output:
!   run       -- fel_run_struct: Beam built (%fbeam, %nslice).
!   err_flag  -- logical: Set True if there is an error. False otherwise.
!-

function fel_line_slippage (run, lambda) result (slip)

type (fel_run_struct), target :: run
type (branch_struct), pointer :: br
real(rp) lambda, slip
integer je

! The line slips one wavelength per undulator period, so this is the headroom a window
! needs at its head for radiation that leaves the bunch behind (fel-physics.md sec-window).

slip = 0
br => run%lat%branch(0)
do je = 1, br%n_ele_track
  if (.not. run%is_fel(je)) cycle
  if (br%ele(je)%value(l_period$) <= 0) cycle
  slip = slip + br%ele(je)%value(l$) / br%ele(je)%value(l_period$)
enddo
slip = slip * lambda

end function fel_line_slippage

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine fel_init_beam (run, err_flag)

type (fel_run_struct), target :: run
logical err_flag

type (lat_struct), pointer :: lat
type (branch_struct), pointer :: branch
type (fel_beam_struct), pointer :: fbeam
type (beam_init_struct), pointer :: beam_init
type (fel_resample_param_struct) resample
character(400) beam_file, write_genesis_dist, write_openpmd_file
character(16) load_mode
character(400) out_root
character(400) field_file(9)
logical quiet_start, shot_noise, gen_test_weights, resample_split_weights
logical split_weights, swap_beam_xy, err
integer beamlet_size, ran_seed, is, ih
real(rp) gamma0, lambda0, window_length, seed_power, seed_waist_size, grid_half_width
integer n_win
integer window_sample, grid_n_pts, npart_gen
real(rp) delgam_gen, tw_beta_x, tw_alpha_x, tw_beta_y, tw_alpha_y
character(*), parameter :: r_name = 'fel_init_beam'

!

err_flag = .false.
lat => run%lat
branch => lat%branch(0)
fbeam => run%fbeam
beam_init => run%beam_init
resample = run%resample
gamma0 = run%gamma0
beam_file = run%bparam%beam_file
write_genesis_dist = run%bparam%write_genesis_dist
write_openpmd_file = run%bparam%write_openpmd_file
load_mode = run%bparam%load_mode
quiet_start = run%bparam%quiet_start
beamlet_size = run%bparam%beamlet_size
shot_noise = run%bparam%shot_noise
gen_test_weights = run%bparam%gen_test_weights
resample_split_weights = run%bparam%resample_split_weights
split_weights = run%bparam%split_weights
swap_beam_xy = run%bparam%swap_beam_xy
ran_seed = run%global%ran_seed
out_root = run%global%out_root
field_file = run%field_file
lambda0 = run%winit%lambda0
window_length = run%slicing%window_length
window_sample = run%slicing%n_wavelength

! A slice count states the window as well as a length does, and Genesis's own decks
! think in one or the other. Both together would have to be reconciled, so they are
! refused instead.

if (run%slicing%n_slice > 0) then
  if (window_length > 0) then
    call out_io (s_error$, r_name, 'SLICING%WINDOW_LENGTH AND SLICING%N_SLICE BOTH STATE THE', &
                                   'WINDOW. GIVE ONE.')
    err_flag = .true.;  return
  endif
  window_length = run%slicing%n_slice * window_sample * lambda0
endif

! A flat current is Genesis's &beam current: the bunch's z structure is then not
! consulted, and with no window it is the steady state of one slice. It and
! beam_init%bunch_charge are two ways to say the same thing, so both is refused.

if (run%slicing%current > 0 .and. beam_init%bunch_charge > 0) then
  call out_io (s_error$, r_name, 'SLICING%CURRENT AND BEAM_INIT%BUNCH_CHARGE BOTH SET THE CHARGE.', &
                                 'GIVE ONE.')
  err_flag = .true.;  return
endif
grid_n_pts = run%winit%grid_n_pts
grid_half_width = run%winit%grid_half_width
seed_power = run%winit%seed_power
seed_waist_size = run%winit%seed_waist_size

! One reference energy, and the lattice is it: gamma0 = e_tot/m_e c^2 from the lattice
! header, never a namelist input. There used to be a namelist gamma0 for Genesis-deck
! symmetry. The first external user fed it a hand-rounded value against a round lattice
! e_tot, and the two disagreed at 1.4e-9. The run died mid-tracking on the seam's
! backstop p0c check with raw numbers: the FEL physics ran on one reference while
! Bmad's momenta were normalized by the other. Two specifications of one truth is the
! defect. The redundant one was removed (parameters live on the lattice).

if ((beam_file == '') .neqv. (field_file(1) == '')) then
  call out_io (s_error$, r_name, 'GIVE BOTH BEAM_FILE AND FIELD_FILE, OR NEITHER (TO GENERATE).')
  err_flag = .true.;  return
endif
if (field_file(1) == '' .and. any(field_file(2:) /= '')) then
  call out_io (s_error$, r_name, 'HARMONIC FIELD FILES NEED THE FUNDAMENTAL IN FIELD_FILE(1).')
  err_flag = .true.;  return
endif
select case (trim(load_mode))
case ('sample', 'keep')
case default
  call out_io (s_error$, r_name, 'LOAD_MODE MUST BE "sample" OR "keep", GOT: ' // trim(load_mode))
  err_flag = .true.;  return
end select
if (beam_file /= '' .and. (beam_init%position_file /= '' .or. resample%use_beam_init)) then
  call out_io (s_error$, r_name, 'A BEAM ALREADY IN SLICES (beam_file) AND A BUNCH TO SLICE', &
                                 '(beam_init%position_file) ARE MUTUALLY EXCLUSIVE.')
  err_flag = .true.;  return
endif

if (beam_file /= '') then
  if (.not. file_is_openpmd(beam_file)) then
    call out_io (s_error$, r_name, 'BEAM FILE IS NOT openPMD: ' // trim(beam_file), &
                 'THIS TRACKER READS openPMD ONLY. CONVERT A GENESIS DUMP FIRST:', &
                 '  lucifer/tests/scripts/convert_genesis.py to-openpmd <in.par.h5> <out.beam.h5>')
    err_flag = .true.;  return
  endif

  ! An openPMD beam file carries the slice partition and not the radiation it was sliced
  ! on, so the deck states the window. Refused rather than defaulted: a wrong wavelength
  ! rescales every phase in the run.
  if (lambda0 <= 0) then
    call out_io (s_error$, r_name, 'READING A BEAM DUMP NEEDS LAMBDA0 > 0.', &
                 'AN openPMD BEAM FILE CARRIES THE SLICE PARTITION AND NOT THE RADIATION', &
                 'IT WAS SLICED ON, SO THE DECK MUST STATE lambda0 AND slicing%n_wavelength.', &
                 'FILE: ' // trim(beam_file))
    err_flag = .true.;  return
  endif

  ! window_length states a window, so the file must match it. Without one the file's own
  ! patch count is the window, which is the usual restart.
  n_win = -1
  if (window_length > 0 .and. window_sample > 0) &
                              n_win = nint(window_length / (window_sample * lambda0))
  call fel_read_openpmd_beam (fbeam, beam_file, gamma0, n_win, lambda0, &
                              window_sample * lambda0, branch%ele(0), err)
  if (err) then
    err_flag = .true.;  return
  endif
elseif (load_mode == 'keep') then
  call keep_initial_state ()
  if (err_flag) return
elseif (beam_init%position_file /= '' .or. resample%use_beam_init) then
  call import_initial_state ()
  if (err_flag) return
else
  call generate_initial_state ()
  if (err_flag) return
endif

if (split_weights) call do_split_weights (fbeam)
run%nslice = size(fbeam%slice)

! Every load path ends here, and the labels are handed out once for all of them. A beam
! read from a dump keeps the labels the file carried, so a restart follows the same
! particles the run before it did.

call fel_assign_ids (fbeam)

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
! Routine to generate the beam from the beam_init description, the sample mode with no
! file (manual sec-loading): matched Gaussian transverse planes on the lattice Twiss, the
! derived per-slice current, beamlets on a uniform ponderomotive phase grid under the
! quiet start and independent particles without it, and the shot noise where asked.
!-

subroutine generate_initial_state ()

type (fel_slice_struct), pointer :: sl
real(rp) p0_mc, ks_l, eg_x, eg_y, u, v, x, xp, y, yp, gam, p_mc, beta, pz, theta0
real(rp) dx_grid, w_part, e0, xg, yg, wsum, w2sum, n_lambda, n_eff, floor_b2, target_b2
real(rp) phi, an, nbl, br, bi
real(rp), allocatable :: theta_work(:), beta_work(:), kick(:), cur_gen(:)
real(rp) nl_min, nl_max, neff_min, neff_max, floor_max, floor_nl
real(rp) spacing_gen, zlen_gen, s_i, q_slice, zbunch_gen, slip_gen
integer ncenter_gen
integer ib, im, ip, mbase, ix, iy, is_g, nslice_gen, ih, nharm, n_clamp, rule, rule_seen, n_group_w
integer, allocatable :: group_w(:)
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
if (npart_gen < 1 .or. beamlet_size < 1 .or. (quiet_start .and. mod(npart_gen, beamlet_size) /= 0)) then
  call out_io (s_error$, r_name, 'BEAM_INIT%N_PARTICLE (MACROPARTICLES PER SLICE HERE) MUST BE A', &
                                 'POSITIVE MULTIPLE OF BEAMLET_SIZE.')
  err_flag = .true.;  return
endif
if (beam_init%a_norm_emit <= 0 .or. beam_init%b_norm_emit <= 0) then
  call out_io (s_error$, r_name, 'BEAM_INIT%A_NORM_EMIT AND %B_NORM_EMIT MUST BE POSITIVE.')
  err_flag = .true.;  return
endif
if (beam_init%bunch_charge <= 0 .and. run%slicing%current <= 0) then
  call out_io (s_error$, r_name, 'THE CHARGE IS UNSTATED: SET BEAM_INIT%BUNCH_CHARGE FOR A BUNCH', &
                                 'WITH A Z STRUCTURE, OR SLICING%CURRENT FOR A FLAT ONE.')
  err_flag = .true.;  return
endif

! The Twiss is the lattice's (one specification of one truth, as with e_tot and the
! import path's init_beam_distribution): read the beginning element, refuse
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
                                 'SLICING%N_WAVELENGTH.')
  err_flag = .true.;  return
endif
if (seed_power > 0 .and. seed_waist_size <= 0) then
  call out_io (s_error$, r_name, 'SEED_WAIST_SIZE MUST BE POSITIVE WHEN SEED_POWER > 0.')
  err_flag = .true.;  return
endif

! The window and the per-slice current derive from the beam_init description (manual
! sec-loading): one bulk bunch, evaluated analytically at the slice centers. The
! default window holds every particle of the described bunch and the line's slippage
! ahead of its head, so radiation that slips forward has somewhere to go before it
! leaves. A Gaussian has no last particle, so its extent is the length beyond which a
! slice would hold less than one electron of the charge, a rule that scales with the
! charge where a fixed number of sigmas does not. A flat bunch takes its grid's extent.
! window_length or n_slice overrides all of it and warns when it clips the bunch.
! sig_z = 0 is the steady state (the whole charge in one slice window) and is refused
! for time-dependent windows.

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

  ! The charge in a slice at s is Q spacing exp(-s^2/2 sig_z^2) / (sqrt(2 pi) sig_z), so
  ! the length that holds every electron is where that falls to one electron's worth.
  ! A description too weak to put one electron in its peak slice keeps the old eight
  ! sigmas, there being no length that satisfies the rule.

  q_slice = beam_init%bunch_charge * spacing_gen / (sqrt(twopi) * beam_init%sig_z * e_charge)
  if (q_slice > 1) then
    zlen_gen = 2 * beam_init%sig_z * sqrt(2 * log(q_slice))
  else
    zlen_gen = 8 * beam_init%sig_z
  endif
else
  zlen_gen = 0                            ! Steady state.
endif

! The line's slippage, one wavelength per undulator period, ahead of the window head,
! which is its high-index end (fel-physics.md sec-window). A stated window says its own
! headroom and gets none added.

zbunch_gen = zlen_gen
slip_gen = 0
if (zlen_gen > 0 .and. window_length <= 0 .and. run%slicing%n_slice <= 0) then
  slip_gen = fel_line_slippage (run, lambda0)
  zlen_gen = zlen_gen + slip_gen
endif

if (window_length > 0) then
  if (zlen_gen == 0 .and. run%slicing%current <= 0 .and. window_length > 1.5_rp * spacing_gen) then
    call out_io (s_error$, r_name, 'SIG_Z = 0 (THE STEADY-STATE DESCRIPTION) IS INVALID FOR A', &
                 'TIME-DEPENDENT WINDOW. GIVE THE BUNCH A LENGTH (SIG_Z, OR A GRID EXTENT),', &
                 'OR STATE SLICING%CURRENT, WHICH IS FLAT AND NEEDS NO Z DESCRIPTION.')
    err_flag = .true.;  return
  endif
  if (window_length < zlen_gen .and. run%slicing%current <= 0) then
    call out_io (s_warn$, r_name, 'window_length = \es10.3\ m CLIPS the described bunch (\es10.3\ m).', &
                 r_array = [window_length, zlen_gen])
  endif
  nslice_gen = max(1, nint(window_length / spacing_gen))
else
  nslice_gen = max(1, nint(zlen_gen / spacing_gen))
endif

if (shot_noise .and. nslice_gen < 2) then
  call out_io (s_error$, r_name, 'SHOTNOISE NEEDS A TIME-DEPENDENT WINDOW, THE SAME RULE AS GENESIS.')
  err_flag = .true.;  return
endif

mbase = npart_gen / beamlet_size
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
fbeam%n_wavelength = window_sample
fbeam%s0 = 0
fbeam%beamlet_size = beamlet_size
fbeam%one4one = .false.

if (allocated(fbeam%slice)) deallocate(fbeam%slice)
allocate (fbeam%slice(nslice_gen), cur_gen(nslice_gen))

! The derived per-slice current: slicing%current where the deck states one, which needs
! no z description and is Genesis's own &beam current. Otherwise from the bunch:
! flat Q*c/extent inside the grid extent. Gaussian
! profile at the slice centers, bunch centered in the window. Steady state = the
! whole charge in the one slice window, I = Q*c/spacing.

if (slip_gen > 0) then
  ncenter_gen = max(1, nint(zbunch_gen / spacing_gen))   ! The derived window: bunch low, headroom at the head.
else
  ncenter_gen = nslice_gen                               ! A stated window centers the bunch it is given.
endif

if (run%slicing%current > 0) then
  cur_gen = run%slicing%current
elseif (flat_z) then
  cur_gen = 0
  do is_g = 1, nslice_gen
    s_i = (is_g - 1) * spacing_gen - (ncenter_gen - 1) * spacing_gen / 2
    if (abs(s_i) <= zbunch_gen / 2) cur_gen(is_g) = beam_init%bunch_charge * c_light / zbunch_gen
  enddo
elseif (zlen_gen > 0) then
  do is_g = 1, nslice_gen
    s_i = (is_g - 1) * spacing_gen - (ncenter_gen - 1) * spacing_gen / 2
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
n_clamp = 0;  rule_seen = 0
nl_min = huge(1.0_rp); nl_max = 0; neff_min = huge(1.0_rp); neff_max = 0; floor_max = 0

do is_g = 1, nslice_gen
  sl => fbeam%slice(is_g)
  call fel_slice_reallocate (sl, npart_gen)
  sl%n = npart_gen
  w_part = cur_gen(is_g) * fbeam%slice_spacing / (c_light * npart_gen)

  ! Quiet start: mbase base samples, each replicated at beamlet_size equally spaced
  ! ponderomotive phases (theta0 spread on a uniform grid within one beamlet spacing),
  ! so bunching harmonics below beamlet_size vanish to roundoff. Weights and coordinates
  ! follow the Genesis chart's own map: z = beta*theta/ks with phi0 = 0,
  ! weight = I*slice_spacing/(c*npart). theta and beta are held in work arrays so noise
  ! can kick the phases before the z conversion.

  ! Without the quiet start there are no beamlets: every particle draws its own five
  ! coordinates and its own phase, which is a load with macroparticle noise at every
  ! harmonic, as a real bunch's would be. mbase then counts particles, not beamlets.

  ip = 0
  do ib = 1, merge(mbase, npart_gen, quiet_start)
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

    if (quiet_start) then
      theta0 = (ib - 0.5_rp) * twopi / (beamlet_size * mbase)
      do im = 0, beamlet_size - 1
        ip = ip + 1
        theta_work(ip) = theta0 + im * twopi / beamlet_size
        beta_work(ip) = beta
        sl%x(ip) = x;   sl%px(ip) = xp
        sl%y(ip) = y;   sl%py(ip) = yp
        sl%pz(ip) = pz
        sl%weight(ip) = w_part
      enddo
    else
      call ran_uniform (u)
      ip = ip + 1
      theta_work(ip) = twopi * u
      beta_work(ip) = beta
      sl%x(ip) = x;   sl%px(ip) = xp
      sl%y(ip) = y;   sl%py(ip) = yp
      sl%pz(ip) = pz
      sl%weight(ip) = w_part
    endif
  enddo

  ! Validation knob: alternate beamlet weights 0.25x/1.75x, charge preserving, uniform
  ! within each beamlet so the quiet cancellation is untouched. Exercises every
  ! weighted-noise path (the asymmetry is strong enough that using a slice-uniform
  ! electron count where the beamlet's charge belongs mis-sets <|b|^2> by 56 percent,
  ! far outside the statistical check). Not a physics input.

  if (gen_test_weights) then
    do ib = 1, mbase
      sl%weight((ib-1)*beamlet_size+1 : ib*beamlet_size) = &
              sl%weight((ib-1)*beamlet_size+1 : ib*beamlet_size) * (1 + 0.75_rp * (-1)**ib)
    enddo
  endif

  ! Noise bookkeeping (fel-physics.md sec-noise): real electrons N_lambda = charge/e, effective
  ! macroparticle number N_eff = (sum w)^2/sum w^2, both per slice.

  wsum = sum(sl%weight(1:npart_gen))
  w2sum = sum(sl%weight(1:npart_gen)**2)
  n_lambda = wsum / e_charge
  n_eff = wsum**2 / w2sum
  nl_min = min(nl_min, n_lambda);  nl_max = max(nl_max, n_lambda)
  neff_min = min(neff_min, n_eff); neff_max = max(neff_max, n_eff)

  ! Zero-current slices (Gaussian tails, outside a flat extent) carry no noise:
  ! Genesis's own zero-current skip, shared with the import.

  if (shot_noise .and. wsum > 0) then

    ! The quiet floor is measured and the noise imposed by the one routine every loader
    ! shares (fel_impose_noise): beamlets here, so fel_fawley_noise runs unchanged and its
    ! draw order, two ran_uniform per harmonic per beamlet, is the one every recorded level
    ! was measured against.

    call fel_impose_noise (sl, theta_work, npart_gen, quiet_start, beamlet_size, grid_half_width, &
                           grid_n_pts, is_g, n_clamp, rule, floor_nl, err_flag)
    if (err_flag) return
    rule_seen = max(rule_seen, rule)
    floor_max = max(floor_max, floor_nl)
  endif

  ! The slice's independent transverse samples: the beamlets, or every particle alone.

  if (quiet_start) then
    sl%m_ind = fel_m_ind (sl, beamlet_groups(npart_gen, beamlet_size), mbase)
  else
    call fel_transverse_groups (sl, group_w, n_group_w)
    sl%m_ind = fel_m_ind (sl, group_w, n_group_w)
  endif

  ! To the stored chart: z = beta*theta/ks with phi0 = 0, beta of the base sample.

  do ip = 1, npart_gen
    sl%z(ip) = beta_work(ip) * theta_work(ip) / ks_l
  enddo
enddo

deallocate (theta_work, beta_work)

if (shot_noise) then
  call fel_report_noise (rule_seen, floor_max, nslice_gen)
  call out_io (s_info$, r_name, '  N_lambda per slice: \es10.3\ to \es10.3\ ', &
               '  N_eff per slice:    \es10.3\ to \es10.3\ ', &
               r_array = [nl_min, nl_max, neff_min, neff_max])
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
! fields in the header table and refuses every other field that is set -- a
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
                                 'QUIET-START GENERATOR (SEE THE HONORED-FIELDS TABLE IN doc/input-reference.md).', &
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
! Routine for the sample mode from a bunch (fel-physics.md sec-import): a bunch_struct,
! read from the openPMD-beamphysics file beam_init%position_file or generated from
! beam_init on the validation route resample%use_beam_init, is resampled into FEL slices
! by the transcribed Genesis importdistribution method (fel_import_mod, where the
! algorithm and its provenance are documented). The RNG-free outputs the exactness checks
! read (the analysis moments and the per-slice current profile) are written at full
! precision.
!-

subroutine import_initial_state ()

type (beam_struct), target :: beam_b
type (bunch_struct), pointer :: bp
real(rp) moments(11)
integer is_g, iu_i
logical err_i
character(400) line
character(*), parameter :: r_name = 'import_initial_state'

!

if (lambda0 <= 0) then
  call out_io (s_error$, r_name, 'IMPORT NEEDS LAMBDA0 > 0.')
  err_flag = .true.;  return
endif
if (window_sample < 1) then
  call out_io (s_error$, r_name, 'SLICING%N_WAVELENGTH MUST BE A POSITIVE INTEGER.')
  err_flag = .true.;  return
endif

! One seed governs the whole import: the bunch generation, the resampler's draws and
! the shot noise. Seeding after generation was the first mutation this path caught in
! development: every run then imports a different bunch. The split-weight and
! thread-determinism checks both fail on what looks like resampler noise.

call ran_seed_put (ran_seed)

call fel_bunch_from_beam_init (beam_b, err_i)
if (err_i) then
  err_flag = .true.;  return
endif

bp => beam_b%bunch(1)

! Check knob: coincident split-weight copies before anything downstream sees the bunch.
! The current profile (weighted sums) and the analysis moments (unweighted, over
! coincident copies) must then be bit-identical to the unsplit run.

if (resample_split_weights) call split_bunch_weights (bp)

if (write_genesis_dist /= '') then
  call fel_write_genesis4_distribution (bp, write_genesis_dist, err_i)
  if (err_i) then
    err_flag = .true.;  return
  endif
  call out_io (s_info$, r_name, 'Wrote Genesis distribution file: ' // trim(write_genesis_dist))
endif

if (write_openpmd_file /= '') then
  call hdf5_write_beam (write_openpmd_file, beam_b%bunch(1:1), .false., err_i, lat)
  if (err_i) then
    err_flag = .true.;  return
  endif
  call out_io (s_info$, r_name, 'Wrote openPMD-beamphysics file: ' // trim(write_openpmd_file))
endif

! On the user's path the per-slice count after copies is beam_init%n_particle and the
! copies are beamlet_size, the same two numbers the generator reads. The resample block's
! own counts are for its validation route, resample%use_beam_init, where the bunch's
! particle count is beam_init%n_particle.

if (.not. resample%use_beam_init) then
  if (beam_init%n_particle < 1) then
    call out_io (s_error$, r_name, 'BEAM_INIT%N_PARTICLE, THE MACROPARTICLES PER SLICE AFTER COPIES,', &
                                   'MUST BE POSITIVE.')
    err_flag = .true.;  return
  endif
  resample%n_particle_per_slice = beam_init%n_particle
  resample%beamlet_size = beamlet_size
endif

call fel_import_bunch (bp, gamma0, lambda0, window_sample * lambda0, resample, quiet_start, shot_noise, &
                       grid_half_width, grid_n_pts, fbeam, err_i, moments)
if (err_i) then
  err_flag = .true.;  return
endif
call out_io (s_info$, r_name, 'Imported into \i0\ slices of \i0\ particles.', &
             i_array = [size(fbeam%slice), resample%n_particle_per_slice])

! The RNG-free instruments the exactness checks read (the analysis moments and the
! per-slice current profile) go to a file, not stdout: stdout is for humans and
! nslice current lines are not (doc/user-guide.md). Full precision, one row per slice.
! Written here, at import time, because load_only stops before tracking.

open (newunit = iu_i, file = trim(out_root) // '.import.txt', action = 'write')
write (iu_i, '(a)') '# The distribution import, at full precision. Machine-readable; stdout is not.'
write (iu_i, '(a, i0, a, i0)') '# nslice = ', size(fbeam%slice), '   npart_per_slice = ', resample%n_particle_per_slice
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
! Subroutine fel_bunch_from_beam_init (beam_b, err_i)
!
! Routine to make the one bunch every slicing path starts from. Bmad's init_beam_distribution
! reads beam_init%position_file when it is set, openPMD or its own ASCII, and generates from
! the beam_init description otherwise, so the two sources meet here and the loaders never
! ask which one they have.
!-

subroutine fel_bunch_from_beam_init (beam_b, err_i)

type (beam_struct), target :: beam_b
logical err_i
integer n_save
character(*), parameter :: r_name = 'fel_bunch_from_beam_init'

!

err_i = .false.
if (beam_init%position_file == '' .and. beam_init%n_particle < 1) then
  call out_io (s_error$, r_name, 'BEAM_INIT%N_PARTICLE MUST BE POSITIVE.')
  err_i = .true.;  return
endif

! Bmad's file read keeps only the first beam_init%n_particle particles when it is set,
! and here the count means the macroparticles per slice, so the read sees zero and
! takes the whole file.

n_save = beam_init%n_particle
if (beam_init%position_file /= '') beam_init%n_particle = 0
beam_init%n_bunch = 1
call init_beam_distribution (branch%ele(0), lat%param, beam_init, beam_b, err_i)
beam_init%n_particle = n_save
if (err_i) return

! A file Bmad reads as several bunches is this tracker's own beam dump, one particle
! patch per slice with each slice's time counted from its own start, so read as a bunch
! its slices would collapse onto one another. The dump path is beam_file.

if (size(beam_b%bunch) > 1) then
  call out_io (s_error$, r_name, 'THE FILE HOLDS \i0\ PARTICLE PATCHES, WHICH IS A BEAM DUMP IN SLICES', &
                                 'AND NOT A BUNCH. LOAD IT WITH beam_file, NOT beam_init%position_file.', &
                                 'FILE: ' // trim(beam_init%position_file), i_array = [size(beam_b%bunch)])
  err_i = .true.;  return
endif

if (beam_init%position_file /= '') then
  call out_io (s_info$, r_name, 'Read \i0\ particles from: ' // trim(beam_init%position_file), &
               i_array = [size(beam_b%bunch(1)%particle)])
else
  call out_io (s_info$, r_name, 'Generated \i0\ particles from beam_init.', &
               i_array = [size(beam_b%bunch(1)%particle)])
endif

end subroutine fel_bunch_from_beam_init

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine keep_initial_state ()
!
! Routine for load_mode = "keep": every live particle of the bunch stays, binned by its
! arrival time into the slices, and the slice's charge is what fell into it. With the
! quiet start each particle becomes a beamlet of beamlet_size copies at weight/beamlet_size
! that share its five other coordinates, phases spread over 2 pi about its own, so the load
! is quiet at every harmonic below beamlet_size and the bunch's moments are the particles'.
! Without it each particle keeps its phase and its weight. No random number is drawn
! before the shot noise, so the same bunch loads the same way every time.
!
! The window follows the bunch, ceiling(extent / spacing) slices from the earliest
! particle, or slicing%n_slice or slicing%window_length where the deck states one, the
! bunch centered and the particles outside counted in a warning.
!-

subroutine keep_initial_state ()

type (beam_struct), target :: beam_b
type (bunch_struct), pointer :: bp
type (coord_struct), pointer :: cp
type (fel_slice_struct), pointer :: sl
real(rp), allocatable :: s_k(:)
real(rp) p0_mc, p_mc, gam, beta_p, ks_l, spacing_k, smin, ttotal, offset, theta0, w_copy
real(rp) floor_max, wsum
integer, allocatable :: isl(:), nfill(:), group_w(:)
integer nalive, nslice_k, ncopy, ip, i, is_g, im, n_out, n_clamp, rule_seen, n_group_w
logical err_i
character(*), parameter :: r_name = 'keep_initial_state'

!

if (lambda0 <= 0) then
  call out_io (s_error$, r_name, 'LOADING NEEDS LAMBDA0 > 0.')
  err_flag = .true.;  return
endif
if (window_sample < 1) then
  call out_io (s_error$, r_name, 'SLICING%N_WAVELENGTH MUST BE A POSITIVE INTEGER.')
  err_flag = .true.;  return
endif
if (beamlet_size < 1) then
  call out_io (s_error$, r_name, 'BEAMLET_SIZE MUST BE POSITIVE.')
  err_flag = .true.;  return
endif

! One seed for the bunch and the noise, as on the sample path.

call ran_seed_put (ran_seed)

call fel_bunch_from_beam_init (beam_b, err_i)
if (err_i) then
  err_flag = .true.;  return
endif
bp => beam_b%bunch(1)

if (resample_split_weights) call split_bunch_weights (bp)

if (write_openpmd_file /= '') then
  call hdf5_write_beam (write_openpmd_file, beam_b%bunch(1:1), .false., err_i, lat)
  if (err_i) then
    err_flag = .true.;  return
  endif
  call out_io (s_info$, r_name, 'Wrote openPMD-beamphysics file: ' // trim(write_openpmd_file))
endif

! Arrival time as a length, tau = -z/beta, from the earliest particle, the import's chart.

nalive = count(bp%particle%state == alive$)
if (nalive < 1) then
  call out_io (s_error$, r_name, 'BUNCH HAS NO LIVE PARTICLES.')
  err_flag = .true.;  return
endif
if (sum(bp%particle%charge, mask = bp%particle%state == alive$) <= 0) then
  call out_io (s_error$, r_name, 'BUNCH HAS ZERO TOTAL CHARGE; NOTHING WOULD LASE.', &
    'AN openPMD FILE WITHOUT CHARGE DATA, OR AN UNSET beam_init%bunch_charge, LOADS DARK.')
  err_flag = .true.;  return
endif

allocate (s_k(size(bp%particle)), isl(size(bp%particle)))
s_k = 0
do ip = 1, size(bp%particle)
  cp => bp%particle(ip)
  if (cp%state /= alive$) cycle
  p_mc = (1 + cp%vec(6)) * cp%p0c / m_electron
  gam = sqrt(p_mc**2 + 1)
  s_k(ip) = -cp%vec(5) * gam / p_mc
enddo
smin = minval(s_k, mask = bp%particle%state == alive$)
ttotal = maxval(s_k, mask = bp%particle%state == alive$) - smin

spacing_k = window_sample * lambda0
if (run%slicing%n_slice > 0) then
  nslice_k = run%slicing%n_slice
elseif (window_length > 0) then
  nslice_k = max(1, nint(window_length / spacing_k))
else
  ! The bunch's own extent, and no slippage added. A loaded window is data: it came from
  ! the run or the code that wrote the dump, and widening it would put this program's
  ! beam in a different window from the one the comparison starts in. A continuation that
  ! wants headroom states it. The described bunch, which this program sizes itself, does
  ! get the slippage.
  nslice_k = max(1, ceiling(ttotal / spacing_k - 1e-9_rp))
endif
offset = (nslice_k * spacing_k - ttotal) / 2       ! The bunch centered in a stated window.
if (run%slicing%n_slice <= 0 .and. window_length <= 0) offset = 0

if (shot_noise .and. nslice_k < 2) then
  call out_io (s_error$, r_name, 'SHOT_NOISE NEEDS A TIME-DEPENDENT WINDOW, THE SAME RULE AS GENESIS.')
  err_flag = .true.;  return
endif

! Bin. A particle outside a stated window is dropped and counted.

allocate (nfill(nslice_k))
nfill = 0;  n_out = 0
do ip = 1, size(bp%particle)
  isl(ip) = 0
  if (bp%particle(ip)%state /= alive$) cycle
  is_g = floor((s_k(ip) - smin + offset) / spacing_k) + 1
  if (is_g < 1 .or. is_g > nslice_k) then
    if (is_g == nslice_k + 1 .and. s_k(ip) - smin + offset <= nslice_k * spacing_k) then
      is_g = nslice_k                                  ! The last particle sits on the edge.
    else
      n_out = n_out + 1;  cycle
    endif
  endif
  isl(ip) = is_g
  nfill(is_g) = nfill(is_g) + 1
enddo
if (n_out > 0) then
  call out_io (s_warn$, r_name, '\i0\ particles fall outside the stated window and are dropped.', &
               i_array = [n_out])
endif

! The beam container and the slices, each sized by what fell into it.

ncopy = merge(beamlet_size, 1, quiet_start)
p0_mc = sqrt(gamma0**2 - 1)
ks_l = twopi / lambda0
fbeam%p0c = p0_mc * m_electron
fbeam%phi0 = 0
fbeam%wavelength = lambda0
fbeam%slice_spacing = spacing_k
fbeam%n_wavelength = window_sample
fbeam%s0 = 0
fbeam%beamlet_size = ncopy
fbeam%one4one = .false.
if (allocated(fbeam%slice)) deallocate (fbeam%slice)
allocate (fbeam%slice(nslice_k))
do is_g = 1, nslice_k
  sl => fbeam%slice(is_g)
  call fel_slice_reallocate (sl, max(1, nfill(is_g) * ncopy))
  sl%n = 0
enddo

! Fill. theta0 is the particle's own phase inside its slice, and the copies stand at
! theta0 + 2 pi m / beamlet_size, so z = beta theta / ks is the stored chart's. The
! transverse momenta are exact, Px / (m c) over p0_mc, where Bmad's vec(2) is Px / P0:
! the resampler carries Genesis's gamma x' instead, which is transcription fidelity.

do ip = 1, size(bp%particle)
  if (isl(ip) == 0) cycle
  cp => bp%particle(ip)
  sl => fbeam%slice(isl(ip))
  p_mc = (1 + cp%vec(6)) * cp%p0c / m_electron
  gam = sqrt(p_mc**2 + 1)
  beta_p = p_mc / gam
  theta0 = modulo(ks_l * (s_k(ip) - smin + offset - (isl(ip) - 1) * spacing_k), twopi)
  w_copy = cp%charge / ncopy
  do im = 0, ncopy - 1
    i = sl%n + 1
    sl%n = i
    sl%x(i) = cp%vec(1);  sl%px(i) = cp%vec(2) * cp%p0c / (m_electron * p0_mc)
    sl%y(i) = cp%vec(3);  sl%py(i) = cp%vec(4) * cp%p0c / (m_electron * p0_mc)
    sl%pz(i) = (p_mc - p0_mc) / p0_mc
    sl%z(i) = beta_p * (theta0 + im * twopi / ncopy) / ks_l
    sl%weight(i) = w_copy
  enddo
enddo

! Noise per slice on the groups the load has, then the independent-sample count. The
! beamlet_size the noise routine gets is the deck's, since it sets the harmonics resolved
! and the floor swept, whether or not this load made beamlets of that size.

n_clamp = 0;  rule_seen = 0;  floor_max = 0
do is_g = 1, nslice_k
  sl => fbeam%slice(is_g)
  if (sl%n == 0) then
    sl%m_ind = 1
    cycle
  endif
  wsum = sum(sl%weight(1:sl%n))
  if (shot_noise .and. wsum > 0) then
    call keep_impose_on_slice (sl, is_g, p0_mc, ks_l, beamlet_size, n_clamp, rule_seen, floor_max)
    if (err_flag) return
  endif
  if (quiet_start) then
    sl%m_ind = fel_m_ind (sl, beamlet_groups(sl%n, ncopy), sl%n / ncopy)
  else
    call fel_transverse_groups (sl, group_w, n_group_w)
    sl%m_ind = fel_m_ind (sl, group_w, n_group_w)
  endif
enddo

if (shot_noise) then
  call fel_report_noise (rule_seen, floor_max, nslice_k)
  if (n_clamp > 0) then
    call out_io (s_warn$, r_name, '\i0\ noise groups had fewer than one real electron (nbl clamped to 1);', &
                 'the noise level there is not physical.', i_array = [n_clamp])
  endif
endif

call out_io (s_info$, r_name, 'Kept \i0\ particles as \i0\ macroparticles in \i0\ slices.', &
             i_array = [nalive - n_out, sum(nfill) * ncopy, nslice_k])

call write_load_record (nslice_k)

end subroutine keep_initial_state

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine keep_impose_on_slice (sl, is_g, p0_mc, ks_l, nbins, n_clamp, rule_seen, floor_max)
!
! Routine to impose the shot noise on one kept slice. The phases come back out of the
! stored z chart, the noise is imposed on the groups the load has, and they go back in
! with the same beta per particle. The chart is z = beta theta / ks, so the round trip is
! exact to roundoff.
!-

subroutine keep_impose_on_slice (sl, is_g, p0_mc, ks_l, nbins, n_clamp, rule_seen, floor_max)

type (fel_slice_struct) sl
real(rp) p0_mc, ks_l, floor_max
integer is_g, nbins, n_clamp, rule_seen

real(rp), allocatable :: theta_k(:), beta_k(:)
real(rp) floor_nl, pm
integer k, rule

!

allocate (theta_k(sl%n), beta_k(sl%n))
do k = 1, sl%n
  pm = (1 + sl%pz(k)) * p0_mc
  beta_k(k) = pm / sqrt(pm**2 + 1)
  theta_k(k) = sl%z(k) * ks_l / beta_k(k)
enddo
call fel_impose_noise (sl, theta_k, sl%n, quiet_start, nbins, grid_half_width, grid_n_pts, &
                       is_g, n_clamp, rule, floor_nl, err_flag)
if (err_flag) return
rule_seen = max(rule_seen, rule)
floor_max = max(floor_max, floor_nl)
do k = 1, sl%n
  sl%z(k) = beta_k(k) * theta_k(k) / ks_l
enddo

end subroutine keep_impose_on_slice

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine split_bunch_weights (bp)
!
! Check knob: every particle of the bunch becomes two coincident copies carrying a third
! and two thirds of its charge, before anything downstream sees it. Every weighted sum,
! and in keep mode every moment, must then be bit-identical to the unsplit run.
!-

subroutine split_bunch_weights (bp)

type (bunch_struct), pointer :: bp
integer n0, ip_g

!

n0 = size(bp%particle)
call reallocate_bunch (bp, 2*n0, save = .true.)
do ip_g = 1, n0
  bp%particle(n0+ip_g) = bp%particle(ip_g)
  bp%particle(n0+ip_g)%charge = 2 * bp%particle(ip_g)%charge / 3
  bp%particle(ip_g)%charge = bp%particle(ip_g)%charge / 3
enddo

end subroutine split_bunch_weights

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine write_load_record (nslice_k)
!
! The RNG-free instruments the exactness checks read: the per-slice current and the
! slice moments of what was loaded, at full precision, one row per slice, in the same
! file the sample path writes. Written at load time because load_only stops before tracking.
!-

subroutine write_load_record (nslice_k)

integer nslice_k
type (fel_slice_struct), pointer :: sl
real(rp) w, m(5)
integer iu_i, is_g, k

!

open (newunit = iu_i, file = trim(out_root) // '.import.txt', action = 'write')
write (iu_i, '(a)') '# The load, at full precision. Machine-readable; stdout is not.'
write (iu_i, '(a, i0)') '# nslice = ', nslice_k
write (iu_i, '(a)') '#  slice   current [A]   n   <x> <px> <y> <py> <pz> (charge weighted)'
do is_g = 1, nslice_k
  sl => fbeam%slice(is_g)
  w = sum(sl%weight(1:sl%n))
  m = 0
  if (w > 0) then
    do k = 1, sl%n
      m = m + sl%weight(k) * [sl%x(k), sl%px(k), sl%y(k), sl%py(k), sl%pz(k)]
    enddo
    m = m / w
  endif
  write (iu_i, '(a, i0, a, es24.15e3, a, i0, 5es24.15e3)') 'slice ', is_g, ' ', &
        c_light * w / fbeam%slice_spacing, ' ', sl%n, m
enddo
close (iu_i)
call out_io (s_info$, r_name, 'Wrote ' // trim(out_root) // '.import.txt (the current profile and slice moments).')

end subroutine write_load_record

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
! Routine to build the field set. Read the fundamental (field_file(1): an openPMD
! EXT_Wavefront) or generate the Gaussian seed
! (wavefront_init: seed_power = 0 is a dark start). Initialize the harmonic entries
! on the fundamental's grid. Fill any from openPMD imports matched by photon energy.
! Check the beam/field window consistency. Needs the beam (fel_init_beam first).
! Errors return through err_flag, and nothing here stops.
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
  if (.not. file_is_openpmd(field_file(1))) then
    call out_io (s_error$, r_name, 'FIELD FILE IS NOT openPMD: ' // trim(field_file(1)), &
                 'THIS TRACKER READS openPMD ONLY. CONVERT A GENESIS DUMP FIRST:', &
                 '  lucifer/tests/scripts/convert_genesis.py to-openpmd <in.fld.h5> <out.wf.h5>')
    err_flag = .true.;  return
  endif
  call fel_read_openpmd_into_field (run, field_file(1), 1, err)
  if (err) then
    err_flag = .true.;  return
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
! carries. The fundamental's grid must match, per the same one-window rule the
! fundamental import obeys.

do ih = 2, n_harm
  call wavefront_init (ffield(ih)%wf, size(wf%Ex,1), size(wf%Ex,2), size(wf%Ex,3), &
                       wf%dx, wf%dy, wf%dz, wf%wavelength / ffield(ih)%harm, 'x', wf%ref_position)
enddo

do is = 2, 9
  if (field_file(is) == '') cycle
  if (.not. file_is_openpmd(field_file(is))) then
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
! E = E0*exp(-r^2/w0^2), intensity 1/e^2 radius w0, integrating to seed_power.
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
