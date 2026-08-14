!+
! Program fel_track_test
!
! FEL tracker validated against Genesis 1.3 Version 4: single-slice steady state
! (deliverable 3) and multi-slice time dependence with slippage (deliverable 4). See
! bsim/fel/README.md and the design brief.
!
! The program walks a Bmad lattice and applies the seam of the design (brief section 4.1):
!
!   - FEL segments are real Bmad wiggler/undulator elements with
!     tracking_method = custom (Bmad's semantics for program-supplied tracking, which
!     this driver is; no name matching). Their FEL parameters come from lattice
!     attributes -- aw from b_max/l_period, helicity from field_calc -- with the brief's
!     7.5 assertions enforced by name at setup. They are stepped in the element's own
!     ds_step (Genesis's delz, now a lattice attribute like everything else) with
!     the transcribed Genesis physics (fel_track_mod): transverse push with natural
!     focusing, RK4 ponderomotive advance, source deposition and FFT field solve. Bmad's
!     own tracking is not used inside them; und_transport = "bmad" swaps the transverse
!     maps for a flattened copy of Bmad's periodic-wiggler kernel, a priced model
!     alternative.
!
!   - Every other element: each slice's bunch is converted from the packed FEL arrays to
!     coord_structs and tracked by Bmad (track1_bunch). The packed arrays ARE Bmad
!     coordinates (see fel_beam_mod), so the conversion is a plain copy; the only phase
!     bookkeeping is one advance of the common reference phase phi0 per element. The
!     radiation field is drifted through free space by wavefront_drift.
!
! Time dependence follows from the starting dumps alone: a multi-slice dump makes a
! time-dependent run (slippage active), a single-slice dump the steady state, with no
! separate switch -- the same rule as Genesis, whose imports carry the time window
! (ImportBeam.cpp). Slippage is precomputed as a schedule over the lattice, transcribing
! Lattice::calcSlippage with Genesis's reference gamma:
!
!   - each undulator step slips dz*(1+aw^2)/(2*gamma0^2*lambda) wavelengths;
!   - interlude elements slip zero, but their lengths accumulate, and when an undulator
!     follows, floor(Lz/(2*gamma0^2*lambda)) + 1 wavelengths of autophasing land on the
!     last interlude element before it (Lattice.cpp:173: "auto phasing would always add
!     some slippage" -- the field record shifts an integer number of wavelengths, the
!     particle phase does not move, exactly as in Genesis where the fractional part is
!     commented out);
!   - a trailing interlude section adds the same fixup to the lattice's last element
!     (Lattice.cpp:193).
!
! The schedule is applied through fel_apply_slippage after each step's field solve and
! before its diagnostics -- Gencore's step order (track beam, track field, slippage,
! diagnostics). The field record rotates rather than moves; everything reading it in time
! order (the per-slice diagnostics here, the final dump) goes through fel_field_index.
!
! The starting state is a pair of Genesis dumps (&write of beam and field), so both codes
! track from bitwise-identical initial conditions. Diagnostics matching Genesis's
! definitions are recorded at the same z positions Genesis records them: once at the start
! and once after every integration step, one step per interlude element -- one row per
! slice per record, in time-window order.
!
! Input is a namelist file:
!
!   &fel_track_params
!     lat_file = "aramis.bmad"                 ! Bmad lattice.
!     beam_file = "Aramis-initial.par.h5"      ! Genesis particle dump to start from.
!     field_file = "Aramis-initial.fld.h5"     ! Genesis field dump to start from.
!     out_root = "fel_td"                      ! Prefix for the output files.
!     interlude_model = "bmad"                 ! "bmad" (the seam, default) or "genesis".
!     und_transport = "genesis"                ! In-undulator transverse maps: transcribed
!                                              !   TrackBeam ("genesis", the tier default)
!                                              !   or the flattened Bmad periodic-wiggler
!                                              !   kernel ("bmad"); priced in the README.
!     split_weights = F                        ! Weight-invariance test mode; see below.
!     write_initial = F                        ! Also dump the initial state (Genesis format).
!     migrate = F                              ! Slice migration (see below).
!     migrate_check = F                        ! Verify phase continuity at each migration.
!   &end
!
! migrate = T moves particles between slices when their ponderomotive phase leaves the
! slice window (fel_migrate_slices: the weighted generalization of Genesis's
! one4one-only localSort, brief 6.4), called serially after every element. OFF BY
! DEFAULT, deliberately: the Genesis-comparison tiers run against Genesis WITHOUT
! one4one, which never migrates, so migration would be a physics-model difference
! inside a transcription-level comparison. Particles leaving the window ends are
! dropped with their charge counted and reported per event and in the end-of-run
! summary (Genesis discards them silently). migrate_check = T additionally verifies at
! every migration that the whole-beam bunching is unchanged by the moves -- the z
! adjustment shifts each mover's phase by an exact multiple of 2*pi*sample, so any
! deviation beyond rounding is a bookkeeping bug -- and reports the worst deviation.
!
! Alternatively, leave beam_file and field_file blank and the program generates its own
! steady-state starting condition -- a quiet-start beam and a Gaussian seed field -- from
! namelist parameters, making a self-contained single run with no Genesis anywhere (see
! bsim/fel/examples). Generation parameters, with the &beam / &field names they mirror:
!
!     lambda0 = 1e-10          ! Radiation wavelength [m] (with dumps it comes from the file).
!     gen_npart = 8192         ! Macroparticles per slice; must be divisible by gen_nbins.
!     gen_nbins = 8            ! Beamlet size of the quiet start.
!     gen_current = 3000       ! Slice current [A].
!     gen_delgam = 1.0         ! Gaussian rms energy spread, units of m_e c^2.
!     gen_ex = 4e-7, gen_ey = 4e-7             ! Normalized emittances [m rad].
!     gen_beta_x = 8.5, gen_alpha_x = -0.70    ! Twiss at the entrance.
!     gen_beta_y = 17.4, gen_alpha_y = 1.40
!     gen_power = 5e3          ! Seed power [W]; 0 gives a dark start (pure SASE).
!     gen_waist_size = 30e-6   ! Seed 1/e^2 intensity radius w0 [m], waist at the entrance.
!     gen_ngrid = 255          ! Transverse grid points per side.
!     gen_dgrid = 2e-4         ! Grid half width [m] (Genesis's dgrid).
!     gen_seed = 12345         ! Random seed, so the example is reproducible.
!     gen_slen = 0             ! Time window [m]. 0: one slice, steady state (the default).
!                              !   Positive: nslice = round(slen/(sample*lambda0)) slices.
!     gen_sample = 1           ! Slice spacing / lambda0 (integer, Genesis's sample).
!     gen_shotnoise = F        ! Impose physical shot noise (time-dependent windows only,
!                              !   the same rule as Genesis's dotime condition).
!     gen_test_weights = F     ! Validation knob: alternate beamlet weights 0.25x/1.75x
!                              !   (charge preserving, uniform within each beamlet) to
!                              !   exercise the weighted-noise paths. Not physics input.
!     load_only = F            ! Generate, write <out_root>-initial dumps, exit without
!                              !   tracking. For the shot-noise statistical gate.
!
! Third way in (deliverable 10): import a particle DISTRIBUTION -- an arbitrary bunch,
! resampled into slices by the transcribed Genesis importdistribution method
! (fel_import_mod, where the algorithm and its provenance live). The bunch comes from
! Bmad's beam_init_struct (use_beam_init = T with a beam_init%... block -- Bmad's
! native equivalent of Genesis's &beam description) or from an openPMD-beamphysics
! file (dist_file). npart/nbins/sample/seed and the seed field reuse the gen_
! parameters above; one seed governs generation, resampling and noise. Knobs, named
! after &importdistribution's where one exists:
!
!     imp%slicewidth = 0.01    ! Sampling window / bunch length (Genesis's slicewidth).
!     imp%nslice = 0           ! 0: round(bunch_length/spacing), Genesis's rule.
!                              ! (Genesis's match/center are NOT ported: a Bmad lattice
!                              !   carries its Twiss and beam_init generates matched
!                              !   bunches already -- see fel_import_mod's header.)
!     use_beam_init = F        ! Generate the bunch from the beam_init%... block.
!     dist_file = ""           ! Or read an openPMD-beamphysics file.
!     write_dist_file = ""     ! Write the bunch as a Genesis &importdistribution
!                              !   input (t/p/x/xp/y/yp + charge, t = -tau/c) -- the
!                              !   shared file of the cross-code gates.
!     write_opmd_file = ""     ! Write the bunch as openPMD-beamphysics.
!     imp_split_weights = F    ! Gate knob: coincident w/3 + 2w/3 copies before import.
!
! The quiet start loads gen_npart/gen_nbins base samples of the transverse and energy
! distributions per slice and replicates each at gen_nbins equally spaced ponderomotive
! phases, so every bunching harmonic below gen_nbins is zero to roundoff. With
! gen_shotnoise = T, physical shot noise is imposed on top, Fawley style, transcribed
! from Genesis's ShotNoise::applyShotNoise and GENERALIZED TO WEIGHTS: per beamlet and
! harmonic h = 1..(gen_nbins-1)/2, every particle of the beamlet gets the phase kick
! -a_h*sin(h*theta + phi), phi uniform, a_h = (2/h)*sqrt(-ln(U)/nbl), where nbl is the
! beamlet's REAL electron count -- its charge over e -- rather than Genesis's ne/mpart
! (identical for uniform weights). Kick algebra: a quiet beamlet acquires
! |b(h)| = h*a_h/2, so <|b(h)|^2> per beamlet is 1/nbl, and the charge-weighted slice
! average is sum(W_j^2/nbl_j)/(sum W_j)^2 = e/sum(W) = 1/N_lambda for ANY cross-beamlet
! weight distribution -- physical noise by construction (brief 6.2). Genesis silently
! clamps nbl < 1 (more macroparticles than electrons); this loader warns when it clamps.
!
! The N_eff discipline (brief 6.2): the loader reports per-beam N_lambda and
! N_eff = (sum w)^2/sum w^2 ranges, and REFUSES to impose noise on a slice whose
! pre-noise quiet floor max_h |b(h)|^2 exceeds 1 percent of the target 1/N_lambda --
! an unquiet representation (weights varying within a beamlet, degraded structure)
! cannot carry noise below its own sampling floor, and imposing on top of it would
! produce a silently wrong startup level. For beams this loader generates the floor is
! roundoff; the guard exists for what future resampled input may bring.
!
! interlude_model selects how the field-free elements are handled. "bmad" is the
! deliverable's architecture: track1_bunch for the particles, the exact theta mapping from
! Bmad's z, wavefront_drift for the field. "genesis" instead uses the transcribed Genesis
! interlude step (fel_track_interlude_genesis) everywhere, which prices what the seam
! changes: with it the whole run should agree with Genesis at transcription level, and the
! difference between the two modes is the transport model difference, measured rather than
! argued about. The slippage schedule is identical in both models.
!
! split_weights = T replaces each imported particle by two copies at identical
! coordinates carrying 1/3 and 2/3 of its weight. Every collective observable -- power,
! bunching, the field itself -- must be identical to the unsplit run, because the
! dynamics is per particle and the sources and reductions are linear in the weight. The
! benchmark harness runs this against the unsplit run to test the weighted paths, which
! nothing Genesis produces can test: the Genesis dump format carries no weights, so a
! Genesis comparison only ever sees the uniform case.
!
! Outputs: <out_root>.diag.txt (one row per slice per record: z, slice, field and beam
! diagnostics), <out_root>-final.fld.h5 and <out_root>-final.par.h5 (Genesis-format dumps
! of the end state, for field-by-field comparison; the field dump is unrotated to time
! order first, as writeFieldHDF5 does).
!-

program fel_track_test

use fel_track_mod
use fel_import_mod
use wavefront_hdf5_mod
use beam_mod

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (fel_beam_struct), target :: fbeam
type (wavefront_struct) wf
type (bunch_struct) bunch
type (fel_und_struct) und
type (fel_slip_struct) slip
type (fel_slice_diag_struct) bdiag

real(rp) :: gamma0 = 0                  ! DERIVED from the lattice e_tot.
logical :: split_weights = .false.
logical :: write_initial = .false.
logical :: migrate = .false., migrate_check = .false.
character(400) :: lat_file = '', beam_file = '', field_file = '', out_root = 'fel_track'
character(16) :: interlude_model = 'bmad'
character(16) :: und_transport = 'genesis'     ! In-undulator transverse maps: the
                                               ! transcribed TrackBeam ("genesis", tier
                                               ! default) or the flattened Bmad periodic
                                               ! kernel ("bmad"); the difference is
                                               ! priced, not assumed (brief 10 step 9).

! Collective effects (deliverable 8), Genesis &wake / &efield names with wake_/sc_
! prefixes. All off by default; see the header.
logical :: wake_on = .false.
real(rp) :: wake_loss = 0, wake_radius = 2.5e-3_rp, wake_conductivity = 0, wake_relaxation = 0
real(rp) :: wake_gap = 0, wake_lgap = 1, wake_hrough = 0, wake_lrough = 1
logical :: wake_roundpipe = .true.
character(8) :: wake_material = ''
real(rp) :: sc_rmax = 0
integer :: sc_ngrid = 100, sc_nz = 0, sc_nphi = 0
logical :: sc_longrange = .false.

real(rp) :: lambda0 = 0                  ! Generation parameters; see the header.
real(rp) :: gen_current = 0, gen_delgam = 0, gen_ex = 0, gen_ey = 0
real(rp) :: gen_beta_x = 0, gen_alpha_x = 0, gen_beta_y = 0, gen_alpha_y = 0
real(rp) :: gen_power = 0, gen_waist_size = 0, gen_dgrid = 0
real(rp) :: gen_slen = 0
integer :: gen_npart = 8192, gen_nbins = 8, gen_ngrid = 255, gen_seed = 12345
integer :: gen_sample = 1
logical :: gen_shotnoise = .false., gen_test_weights = .false., load_only = .false.

! Distribution import (deliverable 10): a bunch_struct -- generated natively from
! Bmad's beam_init_struct, or read from an openPMD-beamphysics file -- resampled into
! FEL slices by the transcribed Genesis method (fel_import_mod). The knobs mirror
! &importdistribution's names with an imp_ prefix; npart/nbins/sample/seed and the
! field come from the gen_ parameters (one field generator for both paths).
type (beam_init_struct) :: beam_init     ! Bmad's native bunch description (&beam_init).
type (fel_import_param_struct) :: imp
logical :: use_beam_init = .false.       ! Generate the bunch from beam_init.
character(400) :: dist_file = ''         ! Or read it from an openPMD-beamphysics file.
character(400) :: write_dist_file = ''   ! Write the bunch as a Genesis DISTRIBUTION
                                         ! file (t/p/x/xp/y/yp + charge, t = -tau/c),
                                         ! the shared input of the cross-code gates.
character(400) :: write_opmd_file = ''   ! Write the bunch as openPMD-beamphysics
                                         ! (hdf5_write_beam), the dist_file round trip.
logical :: imp_split_weights = .false.   ! Gate knob: coincident w/3 + 2w/3 copies
                                         ! BEFORE import; RNG-free outputs must not move.

real(rp), allocatable :: ele_slip(:)     ! Slippage applied after each element's last step [wavelengths].
type (fel_und_struct), allocatable :: und_of(:)   ! Per-element FEL parameters, from lattice attributes.
logical, allocatable :: is_fel(:)                 ! Which tracked elements are FEL segments.
real(rp) z_now, ks, qf, und_slip_step, Lz, gamma0_ref
real(rp) charge_dropped_tot, b_dev_max
integer ie, is, istep, n_arg, iu_diag, iu_nml, nslice, prev_ie, n_moved_tot, iu_wake
logical err, timerun

type (fel_collective_struct) coll

character(400) param_file
character(*), parameter :: r_name = 'fel_track_test'

namelist / fel_track_params / lat_file, beam_file, field_file, out_root, &
                           interlude_model, und_transport, &
                           split_weights, write_initial, lambda0, gen_npart, gen_nbins, &
                           gen_current, gen_delgam, gen_ex, gen_ey, gen_beta_x, gen_alpha_x, &
                           gen_beta_y, gen_alpha_y, gen_power, gen_waist_size, gen_ngrid, &
                           gen_dgrid, gen_seed, gen_slen, gen_sample, gen_shotnoise, &
                           gen_test_weights, load_only, migrate, migrate_check, &
                           wake_on, wake_loss, wake_radius, wake_conductivity, wake_relaxation, &
                           wake_roundpipe, wake_material, wake_gap, wake_lgap, wake_hrough, &
                           wake_lrough, sc_rmax, sc_ngrid, sc_nz, sc_nphi, sc_longrange, &
                           beam_init, imp, use_beam_init, dist_file, write_dist_file, &
                           write_opmd_file, imp_split_weights

! Read parameters.

n_arg = command_argument_count()
if (n_arg /= 1) then
  print '(a)', 'Usage: fel_track_test <param_file>'
  stop 1
endif
call get_command_argument (1, param_file)

open (newunit = iu_nml, file = param_file, status = 'old', action = 'read')
read (iu_nml, nml = fel_track_params)
close (iu_nml)

if (und_transport /= 'genesis' .and. und_transport /= 'bmad') then
  print '(a)', 'fel_track_test: und_transport must be "genesis" or "bmad", got: ' // trim(und_transport)
  stop 1
endif

if (interlude_model /= 'bmad' .and. interlude_model /= 'genesis') then
  print '(a)', 'fel_track_test: interlude_model must be "bmad" or "genesis", got: ' // trim(interlude_model)
  stop 1
endif

! Read the lattice and the starting state: a pair of Genesis dumps (the shared-start
! benchmark methodology), or a self-generated steady-state condition when both file
! names are blank.
!
! FEL elements carry tracking_method = custom, and Bmad's bookkeeping (the reference
! time/energy pass inside bmad_parser, any track1 at the seam) resolves custom tracking
! through track1_custom_ptr: point it at the standard periodic-wiggler kernel, so the
! element behaves as the plain Bmad wiggler it is everywhere EXCEPT inside this
! driver's own FEL walk. In particular the reference time acquires the resonant
! undulation delay of brief 7.5 from Bmad's own code, not from anything written here.

! Both hooks are needed: mat6_calc_method resolves to custom too (auto follows the
! tracking method), and make_mat6 calls through a null make_mat6_custom_ptr otherwise.

track1_custom_ptr => fel_ele_as_wiggler
make_mat6_custom_ptr => fel_mat6_as_wiggler

call bmad_parser (lat_file, lat)
branch => lat%branch(0)

gamma0 = branch%ele(0)%value(e_tot$) / m_electron
print '(a, f0.6, a)', 'fel_track_test: gamma0 = ', gamma0, ' (from the lattice e_tot).'

! ONE reference energy, and the lattice is it: gamma0 = e_tot/m_e c^2 from the lattice
! header, never a namelist input. There used to be a namelist gamma0 for Genesis-deck
! symmetry; the first external user fed it a hand-rounded value against a round lattice
! e_tot, the two disagreed at 1.4e-9, and the run died mid-tracking on the seam's
! backstop p0c check with raw numbers -- the FEL physics ran on one reference while
! Bmad's momenta were normalized by the other. Two specifications of one truth is the
! defect; the redundant one was removed (the deliverable-9 rule: parameters live on
! the lattice).

if ((beam_file == '') .neqv. (field_file == '')) then
  print '(a)', 'fel_track_test: give both beam_file and field_file, or neither (to generate).'
  stop 1
endif
if (beam_file /= '' .and. (dist_file /= '' .or. use_beam_init)) then
  print '(a)', 'fel_track_test: dump files and a distribution import are mutually exclusive.'
  stop 1
endif
if (dist_file /= '' .and. use_beam_init) then
  print '(a)', 'fel_track_test: give dist_file or use_beam_init, not both.'
  stop 1
endif

if (beam_file /= '') then
  call fel_read_genesis4_beam (fbeam, beam_file, gamma0, err)
  if (err) stop 1
  call wavefront_read_genesis4 (wf, field_file, err)
  if (err) stop 1
elseif (dist_file /= '' .or. use_beam_init) then
  call import_initial_state ()
else
  call generate_initial_state ()
endif

if (split_weights) call do_split_weights (fbeam)
nslice = size(fbeam%slice)
ks = twopi / wf%wavelength

! The beam and field must describe the same time window: one field slice per beam slice,
! at the same wavelength. Checked, never assumed (FINDINGS.md section 5).

if (size(wf%Ex, 3) /= nslice) then
  print '(2(a, i0))', 'fel_track_test: beam has ', nslice, ' slices but the field has ', size(wf%Ex, 3)
  stop 1
endif
if (abs(wf%wavelength - fbeam%wavelength) > 1e-12_rp * fbeam%wavelength) then
  print '(a, 2es20.12)', 'fel_track_test: beam and field disagree on the wavelength: ', &
                         fbeam%wavelength, wf%wavelength
  stop 1
endif

! Collective effects (deliverable 8): configure, build the wake kernels, hoist the
! convolution once (Genesis's behavior; recomputed at the migration stride when
! migration can change the currents). The per-slice eloss is written to
! <out_root>.wake.txt so the energy-bookkeeping gate can check the applied loss exactly.

coll%efield%on = (sc_nz >= 1 .or. sc_longrange)
coll%efield%rmax = sc_rmax
coll%efield%ngrid = sc_ngrid
coll%efield%nz = sc_nz
coll%efield%nphi = sc_nphi
coll%efield%longrange = sc_longrange

if (coll%efield%on .and. interlude_model == 'bmad') then
  print '(a)', 'fel_track_test: NOTE space charge acts inside undulators and genesis-model'
  print '(a)', '  interludes only; the Bmad seam''s own space charge is deliverable 9''s domain.'
endif

coll%wake%on = wake_on
coll%wake%loss = wake_loss
coll%wake%radius = wake_radius
coll%wake%conductivity = wake_conductivity
coll%wake%relaxation = wake_relaxation
coll%wake%roundpipe = wake_roundpipe
coll%wake%material = wake_material
coll%wake%gap = wake_gap
coll%wake%lgap = wake_lgap
coll%wake%hrough = wake_hrough
coll%wake%lrough = wake_lrough

if (wake_on) then
  call fel_wake_init (coll%wake, nslice, nint(fbeam%slice_spacing / fbeam%wavelength), &
                      fbeam%wavelength, err)
  if (err) stop 1
  call fel_wake_update (coll%wake, fbeam)
  open (newunit = iu_wake, file = trim(out_root) // '.wake.txt', action = 'write')
  call write_wake_block (0.0_rp)
endif

if (write_initial .or. load_only) then
  call wavefront_write_genesis4 (wf, trim(out_root) // '-initial.fld.h5', err, 'x')
  if (err) stop 1
  call fel_write_genesis4_beam (fbeam, trim(out_root) // '-initial.par.h5', err)
  if (err) stop 1
endif

if (load_only) then
  print '(a)', 'fel_track_test: load_only set; initial state written, no tracking.'
  stop 0
endif

! Time dependence follows from the dumps: more than one slice makes a time-dependent run
! with slippage active; one slice is the steady state and fel_apply_slippage is a no-op.

timerun = (nslice > 1)
slip%timerun = timerun
slip%sample = fbeam%slice_spacing / fbeam%wavelength
gamma0_ref = fel_gamma0(fbeam)

! FEL segments are real Bmad wiggler/undulator elements carrying
! tracking_method = custom -- Bmad's own semantics for "the program supplies the
! tracking", which this driver does. Their FEL parameters come from LATTICE ATTRIBUTES,
! not the namelist (deliverable 9): aw (rms, Genesis's convention) derives from b_max
! and l_period through K = c*b_max/(k_u*m_e c^2) -- reference-energy independent -- with
! aw = K for a helical device and K/sqrt(2) for a planar one; helicity from field_calc;
! Genesis's natural-focusing split kx/ky from the helicity defaults (0.5/0.5 helical,
! 0/1 planar; Bmad's kx roll-off attribute is not yet mapped and must be zero). The 7.5
! assertions are enforced here: a wiggler with zero b_max or l_period would silently get
! factor = 0 in Bmad's own kernel (no resonance, no error), and a fieldmap field_calc
! gets osc_amplitude without focusing -- both are refused by name.

call setup_fel_elements ()

! The rest of the schedule: drift autophasing. Interludes accumulate Lz; the last
! interlude before each undulator gets floor(Lz/(2*gamma0^2*lambda)) + 1 wavelengths
! (Lattice.cpp:171-174, guarded there by Lz > 0). The end-of-lattice fixup
! (Lattice.cpp:191-193) is UNGUARDED in Genesis: the last element always gets
! floor(Lz/(2*gamma0^2*lambda)) + 1, which is +1 even with no trailing interlude at all
! ("autophasing is applied in case for [a] second, succeeding run"). Transcribed as is --
! omitting that +1 leaves the field record one rotation short at the very end, found the
! hard way against the single-segment time-dependent run.

allocate (ele_slip(branch%n_ele_track))
ele_slip = 0
Lz = 0
prev_ie = 0

do ie = 1, branch%n_ele_track
  ele => branch%ele(ie)
  if (ele%value(l$) == 0) cycle
  if (is_fel(ie)) then
    if (Lz > 0 .and. prev_ie > 0) then
      ele_slip(prev_ie) = ele_slip(prev_ie) + floor(Lz / (2 * gamma0_ref**2 * wf%wavelength)) + 1
      Lz = 0
    endif
  else
    Lz = Lz + ele%value(l$)
  endif
  prev_ie = ie
enddo
if (prev_ie > 0) then
  ele_slip(prev_ie) = ele_slip(prev_ie) + floor(Lz / (2 * gamma0_ref**2 * wf%wavelength)) + 1
endif

! Diagnostics file, one row per slice per record at Genesis's record positions, slices in
! time-window order.

open (newunit = iu_diag, file = trim(out_root) // '.diag.txt', action = 'write')
write (iu_diag, '(a, i0)') '# nslice = ', nslice
write (iu_diag, '(a)') '#         z            slice        power         on_axis_intensity        bunching        ' // &
      'bunching_phase        mean_gamma          sigma_gamma           sigma_x               sigma_y' // &
      '               current               n_eff'

n_moved_tot = 0
charge_dropped_tot = 0
b_dev_max = 0

z_now = 0
call write_diag_rows()     ! Initial record, matching Genesis's diag before the first step.

! Walk the lattice.

do ie = 1, branch%n_ele_track
  ele => branch%ele(ie)

  ! Zero length elements (Bmad's end marker, for one) get no step and no diagnostic
  ! record: Genesis's unrolled lattice has no counterpart for them.

  if (ele%value(l$) == 0) cycle

  if (is_fel(ie)) then

    ! FEL segment: Genesis's unroll in the element's OWN step -- Bmad's standard
    ! ds_step/num_steps attributes, whose bookkeeper computes exactly Genesis's
    ! num_steps = round(l/ds_step) (attribute_bookkeeper.f90; there is no namelist
    ! step size, the same rule as every other parameter). Equal steps; slippage after
    ! each step's field solve; any end-of-lattice fixup lands on the last step.

    und = und_of(ie)
    und%bmad_transport = (und_transport == 'bmad')
    und_slip_step = (1 + und%aw**2) / (2 * gamma0_ref**2 * wf%wavelength)
    und%nstep = max(1, nint(ele%value(num_steps$)))
    und%dz = ele%value(l$) / und%nstep

    do istep = 1, und%nstep
      call fel_track_und_step (und, fbeam, wf, slip, coll, err)
      if (err) stop 1
      if (istep == und%nstep) then
        call fel_apply_slippage (slip, wf, und%dz * und_slip_step + ele_slip(ie))
      else
        call fel_apply_slippage (slip, wf, und%dz * und_slip_step)
      endif
      z_now = z_now + und%dz
      if (istep == und%nstep) call do_migrate ()
      call write_diag_rows()
    enddo

  elseif (interlude_model == 'bmad') then

    ! The seam: Bmad tracks each slice's bunch (coordinate copies in and out),
    ! wavefront_drift moves the field (every slice, rotation-invariant), and the common
    ! phase phi0 advances by the reference rate with Genesis's drift surrogate
    ! ks/(2*gamma0^2) as the reference wavenumber.
    !
    ! This slice loop is deliberately SERIAL: track1_bunch parallelizes over particles
    ! internally (track1_bunch_hom, "$OMP parallel do if (thread_safe)", on by default
    ! via global_com%mp_threading_is_safe), so the threads are already busy inside each
    ! call, and parallelizing here as well would nest. The FEL step's parallelism over
    ! slices lives in fel_track_mod.

    do is = 1, nslice
      call fel_slice_to_bunch (fbeam, fbeam%slice(is), ele, bunch, err)
      if (err) stop 1
      call track1_bunch (bunch, ele, err)
      if (err) then
        print '(2a)', 'fel_track_test: tracking error in element ', trim(ele%name)
        stop 1
      endif
      call fel_bunch_to_slice (bunch, ele, fbeam%slice(is), err)
      if (err) stop 1
    enddo

    fbeam%phi0 = fbeam%phi0 + ele%value(l$) * &
                    fel_phi0_rate(ks, ks * 0.5_rp / gamma0_ref**2, fel_p0_mc(fbeam))

    call wavefront_drift (wf, ele%value(l$), err)
    if (err) stop 1

    ! The chamber does not end where the undulator does: the wake's energy loss applies
    ! through seam interludes too, as one kick of the element's length (Genesis applies
    ! it every step, and an interlude is one step).

    call fel_wake_apply (coll%wake, fbeam, ele%value(l$))

    call fel_apply_slippage (slip, wf, ele_slip(ie))

    z_now = z_now + ele%value(l$)
    call do_migrate ()
    call write_diag_rows()

  else

    ! Genesis's own interlude model, transcribed, for pricing what the seam changes.

    qf = 0
    if (ele%key == quadrupole$) qf = ele%value(k1$)
    call fel_track_interlude_genesis (qf, ele%value(l$), fbeam, wf, slip, coll, err)
    if (err) stop 1

    call fel_apply_slippage (slip, wf, ele_slip(ie))

    z_now = z_now + ele%value(l$)
    call do_migrate ()
    call write_diag_rows()
  endif
enddo

close (iu_diag)
if (wake_on) close (iu_wake)

if (migrate) then
  print '(a, i0, a, es12.4, a)', 'fel_track_test: migration moved ', n_moved_tot, &
        ' particles; dropped charge ', charge_dropped_tot, ' C off the window ends.'
  if (migrate_check) then
    print '(a, es10.2)', '  worst whole-beam bunching deviation across migrations: ', b_dev_max
  endif
endif

! Final dumps in Genesis format. The field record is unrotated to time order first --
! time window position is holds record slice 1 + mod(is-1+first, nslice) -- which is what
! writeFieldHDF5.cpp:86 does on the fly.

if (slip%first /= 0) then
  wf%Ex = cshift(wf%Ex, shift = slip%first, dim = 3)
  slip%first = 0
endif

call wavefront_write_genesis4 (wf, trim(out_root) // '-final.fld.h5', err, 'x')
if (err) stop 1

call fel_write_genesis4_beam (fbeam, trim(out_root) // '-final.par.h5', err)
if (err) stop 1

print '(a)', 'fel_track_test done.'
print '(a)', '  ' // trim(out_root) // '.diag.txt'
print '(a)', '  ' // trim(out_root) // '-final.fld.h5'
print '(a)', '  ' // trim(out_root) // '-final.par.h5'

call wavefront_fft_free()

!------------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!+
! Generate the starting condition from the gen_* namelist parameters: quiet-start beam
! slices (one, or a time window of them), optional physical shot noise generalized to
! weights, and a Gaussian seed field (or a dark start at gen_power = 0). See the program
! header for the parameter list, the noise algorithm and its provenance, and the N_eff
! guard. The no-noise single-slice path is arithmetic-identical to the deliverable-4
! loader -- same draw order, same operations -- which the bit-identity anchors rely on.
!-

subroutine generate_initial_state ()

type (fel_slice_struct), pointer :: sl
real(rp) p0_mc, ks_l, eg_x, eg_y, u, v, x, xp, y, yp, gam, p_mc, beta, pz, theta0
real(rp) dx_grid, w_part, e0, xg, yg, wsum, w2sum, n_lambda, n_eff, floor_b2, target_b2
real(rp) phi, an, nbl, br, bi
real(rp), allocatable :: theta_work(:), beta_work(:), kick(:)
real(rp) nl_min, nl_max, neff_min, neff_max, floor_max
integer ib, im, ip, mbase, ix, iy, is_g, nslice_gen, ih, nharm, n_clamp
character(*), parameter :: r_name = 'fel_track_test'

!

if (lambda0 <= 0) then
  print '(a)', 'fel_track_test: generation needs lambda0 > 0.'
  stop 1
endif
if (gen_npart < 1 .or. gen_nbins < 1 .or. mod(gen_npart, gen_nbins) /= 0) then
  print '(a)', 'fel_track_test: gen_npart must be a positive multiple of gen_nbins.'
  stop 1
endif
if (gen_current <= 0 .or. gen_ex <= 0 .or. gen_ey <= 0 .or. gen_beta_x <= 0 .or. gen_beta_y <= 0) then
  print '(a)', 'fel_track_test: gen_current, gen_ex, gen_ey, gen_beta_x and gen_beta_y must be positive.'
  stop 1
endif
if (gen_delgam < 0 .or. gen_power < 0 .or. gen_ngrid < 3 .or. gen_dgrid <= 0 .or. gen_sample < 1) then
  print '(a)', 'fel_track_test: check gen_delgam, gen_power, gen_ngrid, gen_dgrid, gen_sample.'
  stop 1
endif
if (gen_power > 0 .and. gen_waist_size <= 0) then
  print '(a)', 'fel_track_test: gen_waist_size must be positive when gen_power > 0.'
  stop 1
endif

! The window: gen_slen <= 0 is the single-slice steady state; otherwise Genesis's count,
! nslice = round(slen/(sample*lambda0)) (GenTime.cpp:70).

if (gen_slen > 0) then
  nslice_gen = nint(gen_slen / (gen_sample * lambda0))
  if (nslice_gen < 1) nslice_gen = 1
else
  nslice_gen = 1
endif

if (gen_shotnoise .and. nslice_gen < 2) then
  print '(a)', 'fel_track_test: gen_shotnoise needs a time-dependent window (gen_slen), the same rule as Genesis.'
  stop 1
endif

mbase = gen_npart / gen_nbins
if (gen_test_weights .and. mod(mbase, 2) /= 0) then
  print '(a)', 'fel_track_test: gen_test_weights needs an even number of beamlets.'
  stop 1
endif

p0_mc = sqrt(gamma0**2 - 1)
fbeam%p0c = p0_mc * m_electron
fbeam%phi0 = 0
fbeam%wavelength = lambda0
fbeam%slice_spacing = gen_sample * lambda0
fbeam%s0 = 0
fbeam%nbins = gen_nbins
fbeam%one4one = .false.

if (allocated(fbeam%slice)) deallocate(fbeam%slice)
allocate (fbeam%slice(nslice_gen))

call ran_seed_put (gen_seed)

ks_l = twopi / lambda0
eg_x = gen_ex / p0_mc                 ! Normalized emittance to geometric.
eg_y = gen_ey / p0_mc
w_part = gen_current * fbeam%slice_spacing / (c_light * gen_npart)

allocate (theta_work(gen_npart), beta_work(gen_npart))
n_clamp = 0
nl_min = huge(1.0_rp); nl_max = 0; neff_min = huge(1.0_rp); neff_max = 0; floor_max = 0

do is_g = 1, nslice_gen
  sl => fbeam%slice(is_g)
  call fel_slice_reallocate (sl, gen_npart)
  sl%n = gen_npart

  ! Quiet start: mbase base samples, each replicated at gen_nbins equally spaced
  ! ponderomotive phases (theta0 spread on a uniform grid within one beamlet spacing),
  ! so bunching harmonics below gen_nbins vanish to roundoff. Weights and coordinates
  ! follow fel_read_genesis4_beam: z = beta*theta/ks with phi0 = 0,
  ! weight = I*slice_spacing/(c*npart). theta and beta are held in work arrays so noise
  ! can kick the phases before the z conversion.

  ip = 0
  do ib = 1, mbase
    call ran_gauss (u);  call ran_gauss (v)
    x  = sqrt(eg_x * gen_beta_x) * u
    xp = sqrt(eg_x / gen_beta_x) * (v - gen_alpha_x * u)
    call ran_gauss (u);  call ran_gauss (v)
    y  = sqrt(eg_y * gen_beta_y) * u
    yp = sqrt(eg_y / gen_beta_y) * (v - gen_alpha_y * u)

    call ran_gauss (u)
    gam = gamma0 + gen_delgam * u
    p_mc = sqrt(gam**2 - 1)
    beta = p_mc / gam
    pz = (p_mc - p0_mc) / p0_mc

    theta0 = (ib - 0.5_rp) * twopi / (gen_nbins * mbase)

    do im = 0, gen_nbins - 1
      ip = ip + 1
      theta_work(ip) = theta0 + im * twopi / gen_nbins
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
  ! far outside the statistical gate); not a physics input.

  if (gen_test_weights) then
    do ib = 1, mbase
      sl%weight((ib-1)*gen_nbins+1 : ib*gen_nbins) = &
              sl%weight((ib-1)*gen_nbins+1 : ib*gen_nbins) * (1 + 0.75_rp * (-1)**ib)
    enddo
  endif

  ! Bookkeeping the brief's 6.2 demands: real electrons N_lambda = charge/e, effective
  ! macroparticle number N_eff = (sum w)^2/sum w^2, both per slice.

  wsum = sum(sl%weight(1:gen_npart))
  w2sum = sum(sl%weight(1:gen_npart)**2)
  n_lambda = wsum / e_charge
  n_eff = wsum**2 / w2sum
  nl_min = min(nl_min, n_lambda);  nl_max = max(nl_max, n_lambda)
  neff_min = min(neff_min, n_eff); neff_max = max(neff_max, n_eff)

  if (gen_shotnoise) then

    ! The N_eff guard: measure the pre-noise quiet floor. A representation whose floor
    ! is not far below the target 1/N_lambda cannot carry physical noise -- imposing on
    ! top would give a silently wrong startup level. The sweep covers EVERY harmonic the
    ! beamlet structure can resolve (1..gen_nbins-1), not just the imposed ones: an
    ! unquiet weight pattern can park its floor on a harmonic the imposition never
    ! touches (an alternating within-beamlet pattern lands exactly on gen_nbins/2, found
    ! by the guard's own mutation test) and still corrupt the dynamics through the
    ! nonlinear phase evolution.

    target_b2 = 1 / n_lambda
    floor_b2 = 0
    do ih = 1, gen_nbins - 1
      br = 0; bi = 0
      do ip = 1, gen_npart
        br = br + sl%weight(ip) * cos(ih * theta_work(ip))
        bi = bi + sl%weight(ip) * sin(ih * theta_work(ip))
      enddo
      floor_b2 = max(floor_b2, (br**2 + bi**2) / wsum**2)
    enddo
    floor_max = max(floor_max, floor_b2 * n_lambda)

    if (floor_b2 > 0.01_rp * target_b2) then
      print '(a, i0, a)',      'fel_track_test: slice ', is_g, ': the quiet-start floor is not far below the'
      print '(a)',             '  physical shot-noise level -- this representation cannot carry the requested noise.'
      print '(a, es10.2, a, es10.2)', '  max_h |b(h)|^2 = ', floor_b2, '  vs target 1/N_lambda = ', target_b2
      print '(a, es10.2, a, es10.2)', '  N_eff = ', n_eff, '  N_lambda = ', n_lambda
      stop 1
    endif

    ! Fawley-style shot noise: fel_fawley_noise (fel_beam_mod), the ShotNoise
    ! transcription generalized to weights, shared with the distribution import so the
    ! two paths stay one implementation. Draw order is unchanged from when this block
    ! lived inline here (two ran_uniform per harmonic per beamlet, Genesis's loops).

    call fel_fawley_noise (theta_work(1:gen_npart), sl%weight(1:gen_npart), gen_npart, gen_nbins, n_clamp)
  endif

  ! To the stored chart: z = beta*theta/ks with phi0 = 0, beta of the base sample.

  do ip = 1, gen_npart
    sl%z(ip) = beta_work(ip) * theta_work(ip) / ks_l
  enddo
enddo

deallocate (theta_work, beta_work)

if (gen_shotnoise) then
  print '(a, i0, a)',        'fel_track_test: shot noise imposed on ', nslice_gen, ' slices.'
  print '(2(a, es10.3))',    '  N_lambda per slice: ', nl_min, ' to ', nl_max
  print '(2(a, es10.3))',    '  N_eff per slice:    ', neff_min, ' to ', neff_max
  print '(a, es10.2)',       '  worst quiet floor, |b|^2 * N_lambda: ', floor_max
  if (n_clamp > 0) then
    print '(a, i0, a)',      '  WARNING: ', n_clamp, ' beamlet draws had fewer than one real electron'
    print '(a)',             '  (nbl clamped to 1, as Genesis does silently). The noise level in those'
    print '(a)',             '  beamlets is not physical; use fewer macroparticles or more charge.'
  endif
endif

call generate_seed_field (nslice_gen)

end subroutine generate_initial_state

!------------------------------------------------------------------------------

subroutine generate_seed_field (nslice_f)

! The field: a Gaussian seed at its waist in every slice, E = E0*exp(-r^2/w0^2),
! intensity 1/e^2 radius w0, integrating to gen_power; gen_power = 0 is a dark start.
! Grid convention matches Genesis's dgrid: ngrid points spanning +-dgrid,
! dx = 2*dgrid/(ngrid-1), center on axis. Shared by the built-in generator and the
! distribution import (both make their own beam, neither brings a field).

integer nslice_f, ix, iy, is_g
real(rp) dx_grid, e0, xg, yg

!

if (gen_ngrid < 3 .or. gen_dgrid <= 0) then
  print '(a)', 'fel_track_test: check gen_ngrid and gen_dgrid.'
  stop 1
endif
if (gen_power > 0 .and. gen_waist_size <= 0) then
  print '(a)', 'fel_track_test: gen_waist_size must be positive when gen_power > 0.'
  stop 1
endif

dx_grid = 2 * gen_dgrid / (gen_ngrid - 1)
call wavefront_init (wf, gen_ngrid, gen_ngrid, nslice_f, dx_grid, dx_grid, &
                     fbeam%slice_spacing, lambda0, 'x', 0.0_rp)

if (gen_power > 0) then
  e0 = sqrt(4 * (mu_0_vac * c_light) * gen_power / (pi * gen_waist_size**2))
  do iy = 1, gen_ngrid
    yg = (iy - 1) * dx_grid - gen_dgrid
    do ix = 1, gen_ngrid
      xg = (ix - 1) * dx_grid - gen_dgrid
      wf%Ex(ix, iy, 1) = e0 * exp(-(xg**2 + yg**2) / gen_waist_size**2)
    enddo
  enddo
  do is_g = 2, nslice_f
    wf%Ex(:, :, is_g) = wf%Ex(:, :, 1)
  enddo
endif

end subroutine generate_seed_field

!------------------------------------------------------------------------------

subroutine import_initial_state ()

! Deliverable 10: a bunch_struct -- generated from Bmad's beam_init_struct (the native
! equivalent of Genesis's &beam description) or read from an openPMD-beamphysics file
! -- resampled into FEL slices by the transcribed Genesis importdistribution method
! (fel_import_mod, where the algorithm and its provenance are documented). The seed
! field comes from the same generator as the built-in loader. The RNG-free outputs the
! exactness gates read -- the analysis moments and the per-slice current profile --
! are printed at full precision.

type (beam_struct), target :: beam_b
type (bunch_struct), pointer :: bp
real(rp) moments(11)
integer is_g, ip_g, n0
logical err_i

!

if (lambda0 <= 0) then
  print '(a)', 'fel_track_test: import needs lambda0 > 0.'
  stop 1
endif
if (gen_sample < 1) then
  print '(a)', 'fel_track_test: gen_sample must be a positive integer (Genesis''s sample).'
  stop 1
endif

! One seed governs the whole import: the bunch generation, the resampler's draws and
! the shot noise. Seeding AFTER generation was the first mutation this path caught in
! development -- every run then imports a different bunch, and the split-weight and
! thread-determinism gates both fail on what looks like resampler noise.

call ran_seed_put (gen_seed)

if (use_beam_init) then
  if (beam_init%n_particle < 1) then
    print '(a)', 'fel_track_test: beam_init%n_particle must be positive.'
    stop 1
  endif
  beam_init%n_bunch = 1
  call init_beam_distribution (branch%ele(0), lat%param, beam_init, beam_b, err_i)
  if (err_i) stop 1
  print '(a, i0, a)', 'fel_track_test: generated ', size(beam_b%bunch(1)%particle), &
                      ' particles from beam_init.'
else
  call hdf5_read_beam (dist_file, beam_b, err_i, branch%ele(0))
  if (err_i) stop 1
  print '(a, i0, a)', 'fel_track_test: read ', size(beam_b%bunch(1)%particle), &
                      ' particles from: ' // trim(dist_file)
endif

bp => beam_b%bunch(1)

! Gate knob: coincident split-weight copies before anything downstream sees the bunch.
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
  if (err_i) stop 1
  print '(a)', 'fel_track_test: wrote Genesis distribution file: ' // trim(write_dist_file)
endif

if (write_opmd_file /= '') then
  call hdf5_write_beam (write_opmd_file, beam_b%bunch(1:1), .false., err_i, lat)
  if (err_i) stop 1
  print '(a)', 'fel_track_test: wrote openPMD-beamphysics file: ' // trim(write_opmd_file)
endif

imp%npart = gen_npart
imp%nbins = gen_nbins
call fel_import_bunch (bp, gamma0, lambda0, gen_sample * lambda0, imp, fbeam, err_i, moments)
if (err_i) stop 1

print '(a, i0, a, i0, a)', 'fel_track_test: imported into ', size(fbeam%slice), &
                           ' slices of ', gen_npart, ' particles.'
print '(a, 11es24.15e3)', 'import moments (gavg xavg pxavg yavg pyavg ex ey bx by ax ay):', moments
do is_g = 1, size(fbeam%slice)
  print '(a, i0, a, es24.15e3)', 'import current ', is_g, ': ', &
        c_light * sum(fbeam%slice(is_g)%weight(1:fbeam%slice(is_g)%n)) / fbeam%slice_spacing
enddo

call generate_seed_field (size(fbeam%slice))

end subroutine import_initial_state

!------------------------------------------------------------------------------

subroutine do_split_weights (beam)

! Replace each particle by two coincident copies with weights w/3 and 2w/3. The order --
! all first copies, then all second copies -- keeps the original particles' storage
! order, which keeps the RK4 arithmetic per copy identical to the unsplit run.

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

subroutine setup_fel_elements ()

! Recognize FEL segments -- wiggler/undulator elements with tracking_method = custom,
! Bmad's own semantics for program-supplied tracking, which this driver is -- and derive
! their FEL parameters from lattice attributes, enforcing the brief's 7.5 assertions.
! The stored k1x/k1y wiggler attributes are deliberately NOT read (7.5: their helical
! sign disagrees with the tracking locals; nothing here cross-uses them).

type (ele_struct), pointer :: w
integer je
real(rp) kw, kk

allocate (is_fel(branch%n_ele_track), und_of(branch%n_ele_track))
is_fel = .false.

do je = 1, branch%n_ele_track
  w => branch%ele(je)
  if (.not. (w%key == wiggler$ .or. w%key == undulator$)) cycle
  if (w%tracking_method /= custom$) cycle

  ! The 7.5 assertions live in fel_assert_wiggler_sane, ONE authority called from the
  ! track1/mat6 hooks (where they fire first, during the parse) and again here. Keeping
  ! a second copy inline was tried and rejected: redundant assertions mask the removal
  ! of either copy, which defeats mutation testing of the refusal gates.

  call fel_assert_wiggler_sane (w)

  is_fel(je) = .true.
  kw = twopi / w%value(l_period$)

  ! aw (rms, Genesis's convention) from the peak field:
  ! K = c*b_max/(k_u * m_e c^2), exactly and independent of the reference energy;
  ! helical aw = K, planar aw = K/sqrt(2). Focusing split: Genesis's defaults by
  ! helicity (LatticeParser.cpp:328-333), scaled by ku^2 as Genesis's unroll does.

  kk = c_light * w%value(b_max$) / (kw * m_electron)

  und_of(je)%ku = kw
  und_of(je)%helical = (w%field_calc == helical_model$)
  if (und_of(je)%helical) then
    und_of(je)%aw = kk
    und_of(je)%kx = 0.5_rp * kw**2
    und_of(je)%ky = 0.5_rp * kw**2
  else
    und_of(je)%aw = kk / sqrt(2.0_rp)
    und_of(je)%kx = 0
    und_of(je)%ky = kw**2
  endif
enddo

if (.not. any(is_fel)) then
  print '(a)', 'fel_track_test: the lattice has no FEL elements (wiggler/undulator with tracking_method = custom).'
  stop 1
endif

! Elements carrying Bmad wakes are refused by name: the seam tracks one FEL slice at a
! time as its own bunch, so an element wake would act WITHIN single slices only -- the
! bunch-scale wake between slices never accumulates -- and Bmad's tracker additionally
! notes every zero-charge filler slice ("Wakes are on but bunch charge is zero!").
! Near-null physics applied silently is worse than a refusal. The FEL's own wake model
! (wake_on, deliverable 8) covers in-undulator wakes; bunch-scale element wakes at the
! seam are future work.

do je = 1, branch%n_ele_track
  if (associated(branch%ele(je)%wake)) then
    print '(2a)', 'fel_track_test: element carries Bmad wakes, which the slice-at-a-time ', &
                  'seam cannot apply meaningfully (they would act within single slices only): ' &
                  // trim(branch%ele(je)%name)
    print '(a)',  '  Remove the wakes, or use the FEL wake model (wake_on) for in-undulator wakes.'
    stop 1
  endif
enddo

end subroutine setup_fel_elements

!------------------------------------------------------------------------------

subroutine do_migrate ()

! Slice migration at the per-element stride, serial, between the parallel regions (the
! thread gate stays untouched). Called AFTER z_now is advanced, so per-event drop
! reports carry the z of the diagnostic record they precede -- the conservation
! timeline reconstructs exactly from the log. With migrate_check, the whole-beam
! weighted phasor S = sum(w e^{i theta}) must satisfy S_before = S_after + S_dropped to
! rounding: every mover's phase shifts by an exact multiple of 2*pi*sample and a drop
! removes exactly its own term, so any deviation beyond rounding is a bookkeeping bug
! (wrong z adjustment, weight not moved), not statistics.

real(rp) chd, sb_re, sb_im, sa_re, sa_im, d_re, d_im, wsum
integer nm

if (.not. migrate) return

if (migrate_check) call whole_beam_phasor (sb_re, sb_im, wsum)

call fel_migrate_slices (fbeam, ks, nm, chd, d_re, d_im, err)
if (err) stop 1

n_moved_tot = n_moved_tot + nm
charge_dropped_tot = charge_dropped_tot + chd

! Migration changes the current profile, which the wake convolution was hoisted on
! (brief 4.3's premise predates migration): recompute at this stride. Every recompute
! appends a z-stamped block to <out_root>.wake.txt, so "the wake followed the currents"
! is a structural fact a gate can parse without reimplementing the convolution.

if (nm > 0 .and. coll%wake%on) then
  call fel_wake_update (coll%wake, fbeam)
  call write_wake_block (z_now)
endif
if (chd > 0) then
  print '(a, es22.14, a, es22.14, a)', 'fel_track_test: migration dropped ', chd, &
                                       ' C off the window ends at z = ', z_now, ' m.'
endif

if (migrate_check .and. (nm > 0 .or. chd > 0)) then
  call whole_beam_phasor (sa_re, sa_im, wsum)
  if (wsum > 0) then
    b_dev_max = max(b_dev_max, sqrt((sb_re - sa_re - d_re)**2 + (sb_im - sa_im - d_im)**2) / wsum)
  endif
endif

end subroutine do_migrate

!------------------------------------------------------------------------------

subroutine whole_beam_phasor (s_re, s_im, wsum)

! The whole-beam weighted phasor sum(w e^{i theta}) and total weight, all slices.
! Exactly conserved across migration (moves shift phases by 2*pi*sample multiples;
! drops are accounted separately), which is what migrate_check verifies.

real(rp) s_re, s_im, wsum, theta, beta, p0_mc, w
integer is, ip

s_re = 0; s_im = 0; wsum = 0
p0_mc = fel_p0_mc(fbeam)

do is = 1, size(fbeam%slice)
  do ip = 1, fbeam%slice(is)%n
    w = fbeam%slice(is)%weight(ip)
    beta = fel_beta_of(p0_mc, fbeam%slice(is)%pz(ip))
    theta = fbeam%phi0 + ks * fbeam%slice(is)%z(ip) / beta
    s_re = s_re + w * cos(theta)
    s_im = s_im + w * sin(theta)
    wsum = wsum + w
  enddo
enddo

end subroutine whole_beam_phasor

!------------------------------------------------------------------------------

subroutine write_wake_block (z)

! One block of per-slice eloss, z-stamped. Written at the hoisted update and at every
! migration-stride recompute; the energy-bookkeeping and stale-wake gates parse these.

real(rp) z
integer is_w

write (iu_wake, '(a, es22.14)') '# z = ', z
do is_w = 1, nslice
  write (iu_wake, '(i8, es24.16)') is_w, coll%wake%eloss(is_w)
enddo

end subroutine write_wake_block

!------------------------------------------------------------------------------

subroutine write_diag_rows ()

! One row per slice, slices in time-window order: beam slice is against field slice
! fel_field_index(slip, is, nslice), the rotation Genesis applies at Field.cpp:329.

real(rp) power, on_axis
integer is

do is = 1, nslice
  call fel_field_diag (wf, fel_field_index(slip, is, nslice), power, on_axis)
  call fel_slice_diag (fbeam, fbeam%slice(is), ks, bdiag)

  write (iu_diag, '(es24.16, i8, 10es24.16)') z_now, is, power, on_axis, bdiag%bunching, &
        bdiag%bunching_phase, bdiag%mean_gamma, bdiag%sigma_gamma, bdiag%sigma_x, bdiag%sigma_y, &
        bdiag%current, bdiag%n_eff
enddo

end subroutine write_diag_rows

end program fel_track_test
