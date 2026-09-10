!+
! Module fel_setup_mod
!
! The lattice walk-in of the FEL tracker: hooks and lattice attributes, the parse,
! element recognition with its refusals, the harmonics/field-set validation, the
! collective configuration, the slippage/autophasing schedule, geometry breaks and
! the diagnostics setup. Split as fel_setup_lattice (needs only the inputs) and
! fel_setup_schedule (needs the built starting state). Library contract: errors
! return through err_flag, and nothing here stops. All terminal output goes through
! out_io.
!-

module fel_setup_mod

use fel_struct
use fel_io_mod

implicit none

! The lattice-attribute registration is process-global and idempotent: registered
! once, reused by every later run in the same process (the re-entrancy contract).

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_setup_lattice (run, err_flag)
!
! Routine to turn the parsed inputs into a recognized lattice: the tracking hooks
! and FEL lattice attributes (registered once per process), the parse itself,
! element recognition and its refusals (setup_fel_elements), the harmonics-request
! validation and the field-set allocation. Errors return through err_flag, and nothing
! here stops (the library contract).
!
! Input:
!   run       -- fel_run_struct: Run state carrying the parsed namelist inputs and the lattice
!                  file name.
!
! Output:
!   run       -- fel_run_struct: Lattice parsed and FEL elements recognized (%lat, %gamma0,
!                  %is_fel, %und_of, %fel_mode, %fel_spp, %fel_ramp, %two_pol), field set
!                  allocated (%ffield, %n_harm, %n_banked, %esc_id).
!   err_flag  -- logical: Set True if there is an error. False otherwise.
!-

subroutine fel_setup_lattice (run, err_flag)

type (fel_run_struct), target :: run
logical err_flag

type (lat_struct), pointer :: lat
type (branch_struct), pointer :: branch
type (fel_field_struct), pointer :: ffield(:)
integer, pointer :: fel_mode(:), fel_spp(:)
real(rp), pointer :: fel_ramp(:)
type (fel_und_struct), pointer :: und_of(:)
logical, pointer :: is_fel(:)
character(400) lat_file
character(16) interlude_model
character(1) seed_polarization
logical migrate, reference_run, err
integer harmonics(9)
integer n_harm, ih
character(*), parameter :: r_name = 'fel_setup_lattice'

!

err_flag = .false.
lat => run%lat
lat_file = run%lat_file
interlude_model = run%global%interlude_model
seed_polarization = run%winit%seed_polarization
migrate = run%global%migrate
reference_run = run%global%reference_run
harmonics = run%winit%harmonics

! (The FEL tracking mode is a lattice attribute. Its refusals live in
! setup_fel_elements, and the unaveraged-vs-collective refusal follows setup, once
! any_unavg is known.)

if (interlude_model /= 'bmad' .and. interlude_model /= 'genesis') then
  call out_io (s_error$, r_name, 'INTERLUDE_MODEL MUST BE "bmad" OR "genesis", GOT: ' // trim(interlude_model))
  err_flag = .true.;  return
endif

if (run%global%transport_model /= 'bmad' .and. run%global%transport_model /= 'genesis') then
  call out_io (s_error$, r_name, 'TRANSPORT_MODEL MUST BE "bmad" OR "genesis", GOT: ' // &
               trim(run%global%transport_model))
  err_flag = .true.;  return
endif

if (run%global%unaveraged_steps_per_period < 10) then
  call out_io (s_error$, r_name, 'UNAVERAGED_STEPS_PER_PERIOD IS BELOW THE CONVERGENCE FLOOR OF 10, GOT: ' // &
               int_str(run%global%unaveraged_steps_per_period))
  err_flag = .true.;  return
endif

if (run%global%unaveraged_ramp_periods < -1) then
  call out_io (s_error$, r_name, 'UNAVERAGED_RAMP_PERIODS MUST BE POSITIVE, 0 (THE DEFAULT OF 2), OR THE', &
               'HARD-EDGE SENTINEL -1, GOT: ' // int_str(run%global%unaveraged_ramp_periods))
  err_flag = .true.;  return
endif

! (bmad_com%radiation_damping_on / %radiation_fluctuations_on come straight from the
! &fel_params namelist: Bmad's own switches, exposed directly as Tao exposes them.)

if (bmad_com%radiation_fluctuations_on .and. migrate) then
  call out_io (s_error$, r_name, 'RADIATION_FLUCTUATIONS DRAWS ONE KICK PER BEAMLET (THE QUIET START', &
                                 'CANCELS PER BEAMLET), AND SLICE MIGRATION SCRAMBLES BEAMLET GROUPING. PICK ONE.')
  err_flag = .true.;  return
endif

! Read the lattice and the starting state: a pair of openPMD dumps (the shared-start
! benchmark methodology, where the reference code's dumps are converted at the harness
! boundary), or a self-generated steady-state condition when both file names are blank.
!
! An FEL element carries fel_method and whatever tracking_method it says, which is
! bmad_standard, so Bmad's own periodic-wiggler kernel gives it its reference time and
! its transport everywhere outside this driver's FEL walk. The custom hooks below are
! for a wiggler whose tracking_method is custom, which names some other program's
! tracking: Bmad calls through both pointers and jumps through a null one if either is
! unset, so both are supplied and both track the element as the wiggler it is.

! Both hooks are needed: mat6_calc_method resolves to custom too (auto follows the
! tracking method), and make_mat6 calls through a null make_mat6_custom_ptr otherwise.

track1_custom_ptr => fel_ele_as_wiggler
make_mat6_custom_ptr => fel_mat6_as_wiggler

! err_flag matters: bmad_parser reports attribute errors (e.g. a wake on an element
! type that cannot carry one) and returns. Without the check the run continues on a
! partial lattice, found when an example's drift wakes silently never attached.
!
! exit_on_error is cleared first, because a fatal parse error otherwise reaches
! err_exit, which bombs for a traceback and then stops. That stop carries status zero,
! so a refused lattice looked like a successful run to anything reading exit codes. With
! the flag cleared, bmad_parser returns and the library contract holds: nothing here
! stops, the caller decides.

global_com%exit_on_error = .false.
call bmad_parser (lat_file, lat, err_flag = err)
if (err) then
  call out_io (s_error$, r_name, 'LATTICE PARSE ERRORS (ABOVE); REFUSING TO RUN ON A PARTIAL LATTICE.')
  err_flag = .true.;  return
endif
branch => lat%branch(0)

run%gamma0 = branch%ele(0)%value(e_tot$) / m_electron
call out_io (s_info$, r_name, 'gamma0 = \f0.6\ (from the lattice e_tot).', r_array = [run%gamma0])

call setup_fel_elements ()   ! Needs only the parsed lattice. Sets run%two_pol for the
                             ! beam/field construction below.
if (err_flag) return
is_fel => run%is_fel
fel_mode => run%fel_mode

call derive_grid ()          ! Needs the FEL elements, and runs before any field is built.
if (err_flag) return

call derive_slicing_and_step ()   ! Needs the FEL elements, and runs before the beam is sliced.
if (err_flag) return

! The field set (fel-physics.md sec-field-set). Validate the harmonics request, then allocate
! the set and point the fundamental aliases at entry 1 before any field construction:
! every construction path below builds the fundamental through wf.

n_harm = count(harmonics /= 0)
run%n_harm = n_harm
if (harmonics(1) /= 1) then
  call out_io (s_error$, r_name, 'HARMONICS(1) MUST BE 1 -- THE FUNDAMENTAL ANCHORS THE PHASE,', &
                                 'THE PHI0 ADVANCE AND THE SLIPPAGE SCHEDULE; HARMONIC FIELDS RIDE ON IT.')
  err_flag = .true.;  return
endif
do ih = 2, 9
  if (harmonics(ih) == 0) then
    if (any(harmonics(ih:) /= 0)) then
      call out_io (s_error$, r_name, 'HARMONICS MUST BE A GAP-FREE INCREASING LIST (0 PADDING AT THE END).')
      err_flag = .true.;  return
    endif
    exit
  endif
  if (harmonics(ih) <= harmonics(ih-1)) then
    call out_io (s_error$, r_name, 'HARMONICS MUST BE STRICTLY INCREASING (NO DUPLICATES).')
    err_flag = .true.;  return
  endif
enddo
! The load has to carry the harmonic it is asked to start. The quiet start leaves no
! bunching below harmonic beamlet_size, and fel_fawley_noise then imposes the physical
! level on harmonics 1 through (beamlet_size - 1)/2. A harmonic field above that starts
! from a load with no noise at its own frequency, so its SASE has nothing to grow from
! and the run would report a number with no startup behind it. The refusal names the
! beamlet size that would carry it, 2h + 1. A seeded harmonic needs no such noise, so
! the check applies where shot noise is asked for.

if (run%bparam%shot_noise .and. maxval(harmonics) > (run%bparam%beamlet_size - 1) / 2) then
  call out_io (s_error$, r_name, 'HARMONIC \i0\ IS ABOVE THE HIGHEST THE LOAD CARRIES SHOT ' // &
               'NOISE AT, WHICH IS (BEAMLET_SIZE - 1)/2 = \i0\ AT BEAMLET_SIZE = \i0\ .', &
               'POSSIBLE SOLUTION: SET BEAM_INIT%BEAMLET_SIZE TO \i0\ OR MORE.', &
               i_array = [maxval(harmonics), (run%bparam%beamlet_size - 1) / 2, &
                          run%bparam%beamlet_size, 2 * maxval(harmonics) + 1])
  err_flag = .true.;  return
endif

if (n_harm > 1 .and. any(fel_mode == unaveraged$ .and. is_fel)) then
  call out_io (s_error$, r_name, 'HARMONIC FIELDS WITH AN UNAVERAGED ELEMENT ARE NOT IMPLEMENTED', &
                                 '(THE UNAVERAGED MODE CARRIES THE FUNDAMENTAL ENVELOPE ONLY; ITS HARMONIC', &
                                 'COUPLINGS ARE VALIDATED THROUGH THE PARTICLE SPECTRA, NOT A CARRIED FIELD).')
  err_flag = .true.;  return
endif
if (n_harm > 1 .and. run%two_pol) then
  call out_io (s_error$, r_name, 'HARMONIC FIELDS WITH TWO LIVE POLARIZATIONS ARE NOT VALIDATED', &
                                 'TOGETHER YET; RUN ONE OR THE OTHER.')
  err_flag = .true.;  return
endif

! The source model (fel-physics.md sec-coherent-source): validated, then stamped onto every
! FEL element. v1 scope refusals: the coherent source carries the
! Fundamental of one polarization, and it cannot live inside the unaveraged mode,
! which is the independent check (variance reduction inside the explicit-everything
! mode). Even harmonics are invalid
! in the method itself -- F(z,0) = 0 -- and odd ones are a named follow-on.

select case (run%global%source_model)
case ('deposit')
case ('coherent')
  if (any(fel_mode == unaveraged$ .and. is_fel)) then
    call out_io (s_error$, r_name, 'SOURCE_MODEL = "coherent" WITH AN UNAVERAGED ELEMENT: THE', &
                                   'UNAVERAGED MODE RESOLVES EVERYTHING EXPLICITLY, SO VARIANCE REDUCTION IS REFUSED.')
    err_flag = .true.;  return
  endif
  if (n_harm > 1) then
    call out_io (s_error$, r_name, 'SOURCE_MODEL = "coherent" WITH HARMONIC FIELDS IS NOT IN V1', &
                                   '(EVEN HARMONICS ARE INVALID IN THE METHOD; ODD ONES ARE A NAMED FOLLOW-ON).')
    err_flag = .true.;  return
  endif
  if (run%two_pol) then
    call out_io (s_error$, r_name, 'SOURCE_MODEL = "coherent" WITH TWO LIVE POLARIZATIONS IS NOT IN V1.')
    err_flag = .true.;  return
  endif
  if (run%winit%seed_power <= 0 .and. run%field_file(1) == '') then

    ! Measured, not assumed (check_coherent's SASE experiment): the coherent source
    ! understates dark-start startup by ~175x on the reference configuration.
    ! Spontaneous, spatially-incoherent emission dominates SASE startup and is
    ! exactly what the coherent model drops. The slice bunch factor's physical noise
    ! (Fawley, <|B|^2> N_lambda = 1) is present but is not the dominant seed.
    ! Refused. Seeded runs are fully supported.

    call out_io (s_error$, r_name, 'SOURCE_MODEL = "coherent" WITH A DARK START IS REFUSED: THE', &
                 'COHERENT SOURCE DROPS THE SPATIALLY-INCOHERENT SPONTANEOUS EMISSION THAT', &
                 'DOMINATES SASE STARTUP (MEASURED ~175X LOW).', &
                 'POSSIBLE SOLUTION: SEED THE FIELD (SEED_POWER OR FIELD_FILE) OR USE SOURCE_MODEL = "deposit".')
    err_flag = .true.;  return
  endif
  where (is_fel) und_of%source_model = fel_source_coherent$
case default
  call out_io (s_error$, r_name, 'SOURCE_MODEL MUST BE "deposit" OR "coherent", GOT: ' // &
               trim(run%global%source_model))
  err_flag = .true.;  return
end select

! The source filter (fel-physics.md sec-source-filter), validated and stamped onto every
! FEL element the same way. It multiplies the transformed source, which the unaveraged
! mode does not build (it deposits from the resolved motion) and which the coherent source
! has already replaced by an analytic Gaussian carrying no wide-angle content. Neither is
! refused. The filter reaches the elements that build a source to filter and does nothing
! where there is none, as transport_model reaches the averaged elements only: it is on by
! default, so a deck that asks for the unaveraged mode has not asked for the filter and
! must not be stopped by it. A mixed line filters its averaged elements and leaves the
! unaveraged one alone, which is what the per-element filter state already does. Both
! cases say so at setup, since a switch that quietly does nothing is worth a line.
! Genesis turns its own filter off for a non-positive cut or width, since those are
! divisors and a width of zero has no sigmoid (initSourceFilter, FieldSolverFFT.cpp:158-170).
! A run that asked for a filter it can have and silently did not get it is worse than one
! that stops, so those still refuse.

! The coherent source is turned off ahead of the block below rather than inside it. The
! block arms the per-element filter state, and the edge is derived further down under the
! same switch, so a run that fell through to the arming would filter with no edge set.

if (run%global%source_filter .and. run%global%source_model == 'coherent') then
  call out_io (s_info$, r_name, 'The source filter is off for this run: the coherent source is ' // &
               'an analytic Gaussian and carries no wide-angle content to remove.')
  run%global%source_filter = .false.
endif

if (run%global%source_filter) then
  if (any(fel_mode == unaveraged$ .and. is_fel)) then
    call out_io (s_info$, r_name, 'The source filter reaches the averaged elements only. The ' // &
                 'unaveraged mode deposits from the resolved motion and builds no source to filter.')
  endif
  if (run%global%source_filter_width <= 0) then
    call out_io (s_error$, r_name, 'SOURCE_FILTER_WIDTH MUST BE POSITIVE: \es12.3\ ', &
                 r_array = [run%global%source_filter_width])
    err_flag = .true.;  return
  endif

  ! The tolerance fixes the amplitude the sigmoid has to pass at the angle it protects, and
  ! both ends of its range are singular in the margin (fel_filter_margin). At zero the
  ! logarithm of the tolerance diverges. At one the margin runs to positive infinity and the
  ! edge collapses onto the axis, so a caller testing only the margin's sign would take it.
  ! Checked here, before anything takes a logarithm, and whether or not the derived path
  ! will use it.

  if (run%global%source_filter_tolerance <= 0 .or. run%global%source_filter_tolerance >= 1) then
    call out_io (s_error$, r_name, 'SOURCE_FILTER_TOLERANCE MUST LIE STRICTLY BETWEEN ' // &
                 '0 AND 1: \es12.3\ ', &
                 'POSSIBLE SOLUTION: LEAVE IT AT 0.75, WHICH PUTS THE EDGE ON THE PROTECTED ANGLE.', &
                 r_array = [run%global%source_filter_tolerance])
    err_flag = .true.;  return
  endif

  ! The width is a fraction of the edge, so the sigmoid on axis is 1/(1 + exp(-1/width)).
  ! Genesis4's own default of 1 puts that at 0.73, which attenuates the coherent source as
  ! much as the wide angles: measured, it halves the physical in-cone startup power and
  ! delays saturation by five to nine metres, with nothing in the output saying so
  ! (doc/startup-noise.md). A filter meant to remove wide angles keeps the axis, so a
  ! transmission under 0.99 is refused with the width that would pass. The
  ! validation-internal cuts lift the guard, since placing the edge where Genesis4 places
  ! it means taking its width too, and that comparison is the one thing the soft sigmoid
  ! is still for.

  if (run%global%source_filter_xcut == 0 .and. run%global%source_filter_ycut == 0 .and. &
      1 / (1 + exp(-1 / run%global%source_filter_width)) < 0.99_rp) then
    call out_io (s_error$, r_name, 'SOURCE_FILTER_WIDTH = \es12.3\ LEAVES THE SIGMOID AT ' // &
                 '\f6.3\ ON AXIS, SO THE FILTER ATTENUATES THE COHERENT SOURCE.', &
                 'POSSIBLE SOLUTION: USE SOURCE_FILTER_WIDTH AT OR BELOW \es9.2\ , WHICH KEEPS 0.99.', &
                 r_array = [run%global%source_filter_width, &
                            1 / (1 + exp(-1 / run%global%source_filter_width)), 0.217_rp])
    err_flag = .true.;  return
  endif

  ! The edge is either an angle, which the run may derive, or Genesis4's grid-relative
  ! cuts, which exist so the transcription check can place the edge where Genesis4 does.
  ! Both together have no meaning, so they are refused rather than ranked.

  if (run%global%source_filter_angle /= 0 .and. &
      (run%global%source_filter_xcut /= 0 .or. run%global%source_filter_ycut /= 0)) then
    call out_io (s_error$, r_name, 'SOURCE_FILTER_ANGLE AND SOURCE_FILTER_XCUT OR _YCUT ' // &
                 'ARE BOTH SET. THE CUTS ARE VALIDATION-INTERNAL AND MOVE WITH THE GRID.', &
                 'POSSIBLE SOLUTION: SET THE ANGLE ALONE, OR LEAVE IT UNSET TO DERIVE IT.')
    err_flag = .true.;  return
  endif
  if (run%global%source_filter_angle < 0) then
    call out_io (s_error$, r_name, 'SOURCE_FILTER_ANGLE MUST BE POSITIVE: \es12.3\ ', &
                 r_array = [run%global%source_filter_angle])
    err_flag = .true.;  return
  endif

  where (is_fel .and. fel_mode == averaged$)
    und_of%filter%on = .true.
    und_of%filter%width = run%global%source_filter_width
    und_of%filter%mutate = run%global%source_filter_mutate
  end where
endif

allocate (run%ffield(n_harm))
ffield => run%ffield
do ih = 1, n_harm
  ffield(ih)%harm = harmonics(ih)
enddo
allocate (run%n_banked(n_harm), run%esc_id(n_harm))
run%n_banked = 0
run%esc_id = 0

!------------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!+
! Subroutine setup_fel_elements ()
!
! Routine to recognize FEL segments and derive their FEL parameters from lattice
! attributes (Bmad's kx roll-off attribute is not yet mapped and must be zero). An FEL
! segment is a wiggler/undulator tracked by an FEL method: Bmad's own
! semantics for program-supplied tracking, which this driver is. The wiggler sanity
! assertions are enforced: a wiggler with zero b_max or l_period would silently get
! factor = 0 in Bmad's own kernel (no resonance, no error), and a fieldmap field_calc
! gets osc_amplitude without focusing. Both are refused. The stored k1x/k1y
! wiggler attributes are deliberately NOT read: their helical sign disagrees with the
! tracking locals. Nothing here cross-uses them.
!-

subroutine setup_fel_elements ()

type (ele_struct), pointer :: w
integer je
real(rp) kw, kk, rv
logical err_a

allocate (run%is_fel(branch%n_ele_track), run%und_of(branch%n_ele_track))
allocate (run%fel_mode(branch%n_ele_track), run%fel_spp(branch%n_ele_track), run%fel_ramp(branch%n_ele_track))
is_fel => run%is_fel;  und_of => run%und_of
fel_mode => run%fel_mode;  fel_spp => run%fel_spp;  fel_ramp => run%fel_ramp
is_fel = .false.
fel_mode = 0;  fel_spp = 0;  fel_ramp = 0

! The mode is the element's fel_method, which Bmad carries as it carries
! space_charge_method: Averaged or Unaveraged names the multiparticle physics, and the
! element's tracking_method still names how one particle crosses it. A wiggler whose
! fel_method is Off is not this program's element and is left to the seam.

do je = 1, branch%n_ele_track
  w => branch%ele(je)
  if (.not. (w%key == wiggler$ .or. w%key == undulator$)) cycle
  select case (w%fel_method)
  case (averaged$, unaveraged$);  fel_mode(je) = w%fel_method
  case default;                   cycle
  end select

  fel_spp(je) = run%global%unaveraged_steps_per_period
  fel_ramp(je) = run%global%unaveraged_ramp_periods
  if (fel_ramp(je) == -1) fel_ramp(je) = 0

  ! The wiggler sanity assertions live in fel_assert_wiggler_sane, one authority. It is
  ! called from the track1/mat6 hooks (where they fire first, during the parse) and
  ! again here. Keeping a second copy inline was tried and rejected: redundant
  ! assertions mask the removal of either copy, which defeats mutation testing of the
  ! refusal checks.

  call fel_assert_wiggler_sane (w)

  is_fel(je) = .true.
  kw = twopi / w%value(l_period$)

  ! aw (rms, Genesis's convention) from the peak field:
  ! K = c*b_max/(k_u * m_e c^2), exactly and independent of the reference energy.
  ! Helical aw = K, planar aw = K/sqrt(2). Focusing split: Genesis's defaults by
  ! helicity, scaled by ku^2 as Genesis's unroll does (fel-physics.md sec-element).

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

  ! Tilt: the wiggle-plane rotation, planar only. A tilted helical is a no-op that
  ! reads as confusion, refused. The transcribed-Genesis maps know no tilt, refused.
  ! The polarization 2-vector on (Ex, Ey): planar (cos t, sin t); helical (1,-i)/sqrt2.

  und_of(je)%tilt = w%value(tilt_tot$)
  if (und_of(je)%tilt /= 0) then
    if (und_of(je)%helical) then
      call out_io (s_error$, r_name, 'TILT ON A HELICAL FEL ELEMENT IS A ROTATION OF A CIRCULARLY', &
                   'SYMMETRIC FIELD -- A NO-OP THAT READS AS A MISTAKE: ' // trim(w%name))
      err_flag = .true.;  return
    endif
    if (run%global%transport_model == 'genesis' .and. fel_mode(je) == averaged$) then
      call out_io (s_error$, r_name, 'THE TRANSCRIBED-GENESIS MAPS KNOW NO TILT, SINCE GENESIS4 HAS', &
                   'NONE. SET TRANSPORT_MODEL = "bmad" TO TILT: ' // trim(w%name))
      err_flag = .true.;  return
    endif
  endif
  und_of(je)%cos_t = cos(und_of(je)%tilt)
  und_of(je)%sin_t = sin(und_of(je)%tilt)
  if (und_of(je)%helical) then
    und_of(je)%pol = [cmplx(1.0_rp, 0.0_rp, rp), cmplx(0.0_rp, -1.0_rp, rp)] / sqrt(2.0_rp)
  else
    und_of(je)%pol = [cmplx(und_of(je)%cos_t, 0.0_rp, rp), cmplx(und_of(je)%sin_t, 0.0_rp, rp)]
  endif
enddo

run%two_pol = seed_polarization == 'y'
do je = 1, branch%n_ele_track
  if (is_fel(je) .and. und_of(je)%sin_t /= 0) run%two_pol = .true.
enddo
if (seed_polarization /= 'x' .and. seed_polarization /= 'y') then
  call out_io (s_error$, r_name, 'SEED_POLARIZATION MUST BE "x" OR "y".')
  err_flag = .true.;  return
endif

if (.not. any(is_fel) .and. .not. reference_run) then
  call out_io (s_error$, r_name, &
               'THE LATTICE HAS NO FEL ELEMENT: NO WIGGLER OR UNDULATOR CARRIES', &
               'TRACKING_METHOD = FEL_AVERAGED OR FEL_UNAVERAGED. NOTE THAT CUSTOM IS NOT', &
               'ONE OF THEM: IT MEANS ANOTHER PROGRAM''S TRACKING AND IS LEFT TO THE SEAM.', &
               'POSSIBLE SOLUTION: SET REFERENCE_RUN = T FOR A DELIBERATE NO-FEL RUN (BMAD TRACKS EVERYTHING).')
  err_flag = .true.;  return
endif

end subroutine setup_fel_elements

!------------------------------------------------------------------------------
!+
! Subroutine derive_grid ()
!
! Routine to derive the transverse grid from the beam where the deck states none
! (doc/startup-noise.md, Recommendations). The rms beam size comes from the emittances
! the deck states and the matched Twiss the lattice states, averaged over the FEL
! elements by length, since that is the beta the mode sees. Cells are then a seventh of
! it and the half width nine of it, and the point count follows from whichever half
! width is in force, rounded up to a power of two. The rounding only ever refines, more
! points over one half width being smaller cells, and it leaves the count the same on
! every backend, so a derived deck runs on the device and on the CPU at one grid. A stated value is never overridden, and the other is derived
! against it, so a deck that fixes the cell size gets the containment it needs and one
! that fixes the point count gets the resolution.
!
! Nothing is derived without an emittance to derive from. A run that loads its beam from
! a dump has none, and the refusals that already ask for the grid still ask.
!-

subroutine derive_grid ()

real(rp) sig, dx_want, half
integer n_pts, n_raw

!

if (run%winit%grid_n_pts > 0 .and. run%winit%grid_half_width > 0) return

sig = described_beam_size ()
if (sig <= 0) return

half = run%winit%grid_half_width
if (half <= 0) then
  half = fel_widths_per_sigma$ * sig
  run%winit%grid_half_width = half
  call out_io (s_info$, r_name, 'Grid half width \es10.3\ m, derived: \f0.1\ rms beam sizes ' // &
               'of \es10.3\ m (doc/startup-noise.md).', &
               r_array = [half, fel_widths_per_sigma$, sig])
endif

if (run%winit%grid_n_pts <= 0) then
  dx_want = sig / fel_cells_per_sigma$
  n_raw = 2 * ceiling(half / dx_want) + 1
  n_pts = 1
  do while (n_pts < n_raw)
    n_pts = 2 * n_pts
  enddo
  run%winit%grid_n_pts = n_pts
  call out_io (s_info$, r_name, 'Grid \i0\ points, derived: cells of \es10.3\ m, a beam size ' // &
               'over \f0.1\ , rounded up from \i0\ to a power of two.', &
               i_array = [n_pts, n_raw], r_array = [2 * half / (n_pts - 1), fel_cells_per_sigma$])
endif

end subroutine derive_grid

!------------------------------------------------------------------------------
! contains
!+
! Function described_beam_size () result (sig)
!
! The rms transverse size of the beam the deck describes, from the emittances it states
! and the matched Twiss the lattice states, averaged over the FEL elements by length
! since that is the beta the mode sees. The quadratic mean of the two planes, one size
! for a round description. Zero where there is nothing to derive from, which a run that
! loads its beam from a dump has, and the refusals that already ask still ask.
!
! Output:
!   sig -- real(rp): The rms size [m], or zero.
!-

function described_beam_size () result (sig)

real(rp) sig, eps_a, eps_b, beta_a, beta_b, wt
integer je

!

sig = 0
if (run%gamma0 <= 0) return
eps_a = run%beam_init%a_norm_emit / run%gamma0
eps_b = run%beam_init%b_norm_emit / run%gamma0
if (eps_a <= 0 .or. eps_b <= 0) return

! A lattice whose Twiss was never propagated leaves these at zero, and the beginning
! element is then the one truth there is.

beta_a = 0;  beta_b = 0;  wt = 0
do je = 1, branch%n_ele_track
  if (.not. is_fel(je)) cycle
  if (branch%ele(je)%a%beta <= 0 .or. branch%ele(je)%b%beta <= 0) cycle
  beta_a = beta_a + branch%ele(je)%value(l$) * branch%ele(je)%a%beta
  beta_b = beta_b + branch%ele(je)%value(l$) * branch%ele(je)%b%beta
  wt = wt + branch%ele(je)%value(l$)
enddo
if (wt > 0) then
  beta_a = beta_a / wt;  beta_b = beta_b / wt
else
  beta_a = branch%ele(0)%a%beta;  beta_b = branch%ele(0)%b%beta
endif
if (beta_a <= 0 .or. beta_b <= 0) return

sig = sqrt(0.5_rp * (eps_a * beta_a + eps_b * beta_b))

end function described_beam_size

!------------------------------------------------------------------------------
! contains
!+
! Function described_pierce () result (rho)
!
! The one-dimensional Pierce parameter of the beam the deck describes, at the first FEL
! element, in Genesis's form: rho^3 = (I/I_A) fc^2 / (8 gamma^3 sigma^2 ku^2), with the
! undulator's own coupling, which already carries aw. This is the same expression the
! source filter's edge uses (fel_filter_angles), computed from the description rather
! than from the built beam, because the slice spacing has to be known before the beam is
! sliced at it.
!
! The peak current is the deck's own: a stated flat current, else the Gaussian peak
! Q c / (sqrt(2 pi) sig_z), else the flat extent, else the steady state's whole charge in
! one slice at the stated spacing.
!
! Output:
!   rho -- real(rp): The Pierce parameter, or zero where there is nothing to derive from.
!-

function described_pierce () result (rho)

type (fel_und_struct), pointer :: und
real(rp) rho, sig, cur, fc, zlen, spacing
integer je, je_first

!

rho = 0
sig = described_beam_size ()
if (sig <= 0) return

je_first = 0
do je = 1, branch%n_ele_track
  if (is_fel(je)) then
    je_first = je
    exit
  endif
enddo
if (je_first == 0) return
und => run%und_of(je_first)
if (und%ku <= 0) return
fc = fel_und_coupling (und, 1)
if (fc == 0) return

! The peak current of the description.

if (run%slicing%current > 0) then
  cur = run%slicing%current
else if (run%beam_init%bunch_charge <= 0) then
  return
else if (run%beam_init%sig_z > 0) then
  cur = run%beam_init%bunch_charge * c_light / (sqrt(twopi) * run%beam_init%sig_z)
else
  zlen = run%beam_init%grid(3)%x_max - run%beam_init%grid(3)%x_min
  if (zlen > 0) then
    cur = run%beam_init%bunch_charge * c_light / zlen
  else
    spacing = max(1, run%slicing%n_wavelength) * run%winit%lambda0
    if (spacing <= 0) return
    cur = run%beam_init%bunch_charge * c_light / spacing
  endif
endif
if (cur <= 0) return

rho = ((cur / (4 * pi * m_electron / (mu_0_vac * c_light))) * fc**2 / &
       (8 * run%gamma0**3 * sig**2 * und%ku**2)) ** (1.0_rp / 3.0_rp)

end function described_pierce

!------------------------------------------------------------------------------
! contains
!+
! Subroutine derive_slicing_and_step ()
!
! Routine to derive the slice spacing and the integration step from the gain where the
! deck states neither (fel-physics.md sec-window and sec-element).
!
! The spacing is the whole number of wavelengths in a quarter of the cooperation length
! lambda/(4 pi rho), which is four slices per cooperation length, the resolution FLASH1's
! pulse energy converged at. It is derived only for a time-dependent description. In the
! steady state the spacing is not a resolution at all: the whole charge sits in the one
! slice, so the spacing sets the current, and a derived one would rewrite the beam.
!
! The step is a twentieth of the one-dimensional gain length lambda_u/(4 pi sqrt(3) rho),
! rounded down to whole periods and at least one. Bmad's bookkeeper gives any wiggler
! with a period a step of l_period/20 when the deck states neither ds_step nor num_steps,
! and fills num_steps from it, which resolves the wiggle the averaged model has already
! averaged over. The derived step replaces it. The unset state is that exact quotient,
! since the bookkeeper leaves nothing else to read, so an element whose step is
! l_period/20 by the deck's own choice is derived over. Any other stated step is kept.
!-

subroutine derive_slicing_and_step ()

type (ele_struct), pointer :: ele
real(rp) rho, l_coop, l_1d, step, l_per
integer je, n_wl, n_step, n_done
logical td

!

rho = described_pierce ()

! The spacing. It is derived only for a time-dependent description this program sizes
! itself: in the steady state the whole charge sits in the one slice, so the spacing
! states the current rather than a resolution, and a run that loads its beam has no
! description to derive from. Both of those keep the one wavelength that was the default
! before anything was derived, and the spacing is never left unset.

if (run%slicing%n_wavelength <= 0) then
  td = (run%beam_init%sig_z > 0 .or. &
        run%beam_init%grid(3)%x_max > run%beam_init%grid(3)%x_min .or. &
        run%slicing%window_length > 0 .or. run%slicing%n_slice > 0)
  if (td .and. rho > 0) then
    l_coop = run%winit%lambda0 / (4 * pi * rho)
    n_wl = max(1, int(0.25_rp * l_coop / run%winit%lambda0))
    run%slicing%n_wavelength = n_wl
    call out_io (s_info$, r_name, 'Slice spacing \i0\ wavelengths, derived: a quarter of the ' // &
                 'cooperation length \es10.3\ m at rho = \es10.3\ .', &
                 i_array = [n_wl], r_array = [l_coop, rho])
  else if (td) then
    run%slicing%n_wavelength = 1
    call out_io (s_info$, r_name, 'Slice spacing 1 wavelength: nothing to derive it from, ' // &
                 'the beam being loaded rather than described.')
  else
    run%slicing%n_wavelength = 1
    call out_io (s_info$, r_name, 'Slice spacing 1 wavelength: the steady state holds the whole ' // &
                 'charge in one slice, so the spacing states the current.')
  endif
endif

if (rho <= 0) return

! The step, per element, where the deck states none.

l_1d = 0
n_done = 0
do je = 1, branch%n_ele_track
  if (.not. is_fel(je)) cycle
  ele => branch%ele(je)
  l_per = ele%value(l_period$)
  if (l_per <= 0) cycle
  if (ele%value(ds_step$) /= l_per / 20) cycle

  l_1d = twopi / (run%und_of(je)%ku * 4 * pi * sqrt(3.0_rp) * rho)
  step = max(1.0_rp, real(int(l_1d / 20 / l_per), rp)) * l_per
  n_step = max(1, nint(ele%value(l$) / step))
  ele%value(num_steps$) = n_step
  ele%value(ds_step$) = ele%value(l$) / n_step
  n_done = n_done + 1
enddo

if (n_done > 0) then
  call out_io (s_info$, r_name, 'Integration step \es10.3\ m on \i0\ element(s), derived: a ' // &
               'twentieth of the gain length \es10.3\ m, in whole periods.', &
               i_array = [n_done], r_array = [step, l_1d])
endif

end subroutine derive_slicing_and_step

end subroutine fel_setup_lattice

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_setup_schedule (run, err_flag)
!
! Routine to build everything that needs both the lattice and the built starting state:
! the collective configuration and wake kernels, the wake-window and unaveraged-collective
! refusals, the slippage/autophasing schedule, the geometry breaks (chicanes), and the
! diagnostics setup (dump locators and the exact record/element-end counts). Errors
! return through err_flag, and nothing here stops.
!
! Input:
!   run       -- fel_run_struct: Run state after fel_setup_lattice, fel_init_beam and
!                  fel_init_wavefront (the schedule needs the built beam and field).
!
! Output:
!   run       -- fel_run_struct: Collective configuration (%coll), slippage/autophasing schedule
!                  (%ele_slip, %fel_zoff, %light_corr, %phase_rate, %gamma0_ref), tracking window
!                  (%i_start, %i_end) and diagnostics (%stats, dump locators, record counts) set.
!   err_flag  -- logical: Set True if there is an error. False otherwise.
!-

subroutine fel_setup_schedule (run, err_flag)

type (fel_run_struct), target :: run
logical err_flag

type (lat_struct), pointer :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele, wake_src
type (fel_field_struct), pointer :: ffield(:)
type (fel_beam_struct), pointer :: fbeam
type (wavefront_struct), pointer :: wf
type (fel_collective_struct), pointer :: coll
type (fel_stats_struct), pointer :: stats
integer, pointer :: fel_mode(:)
real(rp), pointer :: ele_slip(:), fel_zoff(:), light_corr(:)
logical, pointer :: is_fel(:), dump_beam_here(:), dump_field_here(:)
character(400) out_root
character(16) interlude_model
character(60) dump_beam_at(40), dump_field_at(40)
logical wake_on, sc_longrange, timerun, two_pol, err
real(rp) sc_rmax, wake_loss, wake_radius, wake_conductivity, wake_relaxation
real(rp) wake_gap, wake_lgap, wake_hrough, wake_lrough
logical wake_roundpipe
character(8) wake_material
integer sc_ngrid, sc_nz, sc_nphi, n_sc, n_ask
real(rp) gamma0_ref, phase_rate, Lz
integer nslice, n_harm, ie, ih, je, prev_ie, nrec_stats, nend_stats
integer harmonics(9)
character(*), parameter :: r_name = 'fel_setup_schedule'

!

err_flag = .false.
lat => run%lat
branch => lat%branch(0)
ffield => run%ffield
fbeam => run%fbeam
wf => run%ffield(1)%wf
coll => run%coll
stats => run%stats
fel_mode => run%fel_mode
is_fel => run%is_fel
out_root = run%global%out_root
interlude_model = run%global%interlude_model
dump_beam_at = run%global%dump_beam_at
dump_field_at = run%global%dump_field_at
harmonics = run%winit%harmonics
nslice = run%nslice
n_harm = run%n_harm
two_pol = run%two_pol
wake_on = run%chamber_wake%on
wake_loss = run%chamber_wake%loss
wake_radius = run%chamber_wake%radius
wake_conductivity = run%chamber_wake%conductivity
wake_relaxation = run%chamber_wake%relaxation
wake_roundpipe = run%chamber_wake%roundpipe
wake_material = run%chamber_wake%material
wake_gap = run%chamber_wake%gap
wake_lgap = run%chamber_wake%lgap
wake_hrough = run%chamber_wake%hrough
wake_lrough = run%chamber_wake%lrough
sc_rmax = run%space_charge%rmax
sc_ngrid = run%space_charge%ngrid
sc_nz = run%space_charge%nz
sc_nphi = run%space_charge%nphi
sc_longrange = run%space_charge%longrange

! Each collective family names its implementation. One value is accepted today, and an
! unknown one is refused rather than silently treated as the transcribed solver.

if (run%chamber_wake%model /= 'genesis') then
  call out_io (s_error$, r_name, 'CHAMBER_WAKE%MODEL MUST BE "genesis", GOT: ' // trim(run%chamber_wake%model))
  err_flag = .true.;  return
endif

if (run%space_charge%model /= 'genesis') then
  call out_io (s_error$, r_name, 'SPACE_CHARGE%MODEL MUST BE "genesis", GOT: ' // trim(run%space_charge%model))
  err_flag = .true.;  return
endif

! The solver's numbers are run-level. Whether it runs in a given element is that
! element's own space_charge_method, resolved into run%sc_here below: nothing here
! decides, so a stats file cannot record a configuration the run did not use.

coll%efield%rmax = sc_rmax
coll%efield%ngrid = sc_ngrid
coll%efield%nz = sc_nz
coll%efield%nphi = sc_nphi
coll%efield%longrange = sc_longrange
coll%efield%model = run%space_charge%model
coll%efield%active = .false.    ! Set per element by the walk, from run%sc_here.

if (sc_nz < 1 .and. .not. sc_longrange .and. bmad_com%csr_and_space_charge_on) then
  ! Neither term is configured, so slice would compute an exact zero at full cost. That
  ! is a deck saying two things at once, and the useful answer is which one it meant.
  ! Only when the master switch is on: with it off, a lattice that carries slice for its
  ! own reasons is being run without space charge deliberately, and demanding solver
  ! numbers for a solver that will not run would be hostile.
  do je = 1, branch%n_ele_track
    if (branch%ele(je)%space_charge_method /= slice$) cycle
    call out_io (s_error$, r_name, &
                 'SPACE_CHARGE_METHOD = SLICE IS SET ON ELEMENT ' // trim(branch%ele(je)%name) // ',', &
                 'BUT NEITHER SPACE-CHARGE TERM IS CONFIGURED: SET SPACE_CHARGE%NZ >= 1 FOR THE', &
                 'SHORT-RANGE HARMONICS, SPACE_CHARGE%LONGRANGE = T FOR THE LONG-RANGE TERM, OR BOTH.')
    err_flag = .true.;  return
  enddo
endif

! Per element: the attribute says whether, Bmad's master switch says whether at all, and
! an FEL element refuses the methods whose solvers need a field this walk does not have.

allocate (run%sc_here(branch%n_ele_track))
run%sc_here = .false.
n_sc = 0

do je = 1, branch%n_ele_track
  select case (branch%ele(je)%space_charge_method)
  case (off$)
    cycle
  case (slice$)
    run%sc_here(je) = bmad_com%csr_and_space_charge_on
    if (run%sc_here(je)) n_sc = n_sc + 1
  case default
    if (run%is_fel(je)) then
      call out_io (s_error$, r_name, 'SPACE_CHARGE_METHOD ON AN FEL ELEMENT MUST BE OFF OR SLICE.', &
                   'THE SLICE-BINNED SOLVE IS THE ONE THIS WALK HAS, SINCE THE FEL SLICES ARE ITS', &
                   'BINS. GOT ' // trim(space_charge_method_name(branch%ele(je)%space_charge_method)) // &
                   ' ON: ' // trim(branch%ele(je)%name))
      err_flag = .true.;  return
    endif
  end select
enddo

! Loud rather than refusing: the master switch exists to disable space charge across a
! lattice without editing it, so the combination is legitimate and only silence is not.

if (.not. bmad_com%csr_and_space_charge_on) then
  n_ask = count(branch%ele(1:branch%n_ele_track)%space_charge_method == slice$)
  if (n_ask > 0) then
    call out_io (s_info$, r_name, &
                 'Note: \i0\ element(s) set space_charge_method = slice, but ' // &
                 'bmad_com%csr_and_space_charge_on', 'is false, so no space charge runs ' // &
                 'anywhere. That combination is how a lattice runs without it.', &
                 i_array = [n_ask])
  endif
endif

if (n_sc > 0 .and. interlude_model == 'bmad') then
  call out_io (s_info$, r_name, 'Note: inside the Bmad seam, space charge is Bmad''s own, from the', &
               'same element attribute. This solver acts in FEL elements and genesis-model interludes.')
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
  if (err) then
    err_flag = .true.;  return
  endif
  call fel_wake_update (coll%wake, fbeam)
  open (newunit = run%iu_wake, file = trim(out_root) // '.wake.txt', action = 'write')
  call fel_write_wake_block (run, 0.0_rp)

  ! The single-particle kernels, written for cross-validation: the same physical
  ! wake fed to Bmad's z_long machinery must reproduce these kernels'
  ! convolution. Note the s = 0 entries carry the Bane self-slice half factor
  ! (fel_wake_init halves them). A plain W(z) table wants the unhalved value.

  if (run%chamber_wake%write_kernels /= '') then
    block
      integer iu_k, i_k
      open (newunit = iu_k, file = trim(run%chamber_wake%write_kernels), action = 'write')
      write (iu_k, '(a)') '# s [m]   wakeres   wakegeo   wakerou   [eV/(m electron)]; s=0 rows are HALVED (Bane self-slice)'
      do i_k = 1, coll%wake%ns
        write (iu_k, '(4es24.15e3)') coll%wake%ds * (i_k-1), coll%wake%wakeres(i_k), &
                                     coll%wake%wakegeo(i_k), coll%wake%wakerou(i_k)
      enddo
      close (iu_k)
      call out_io (s_info$, r_name, 'Wrote wake kernels: ' // trim(run%chamber_wake%write_kernels))
    end block
  endif
endif

! The source filter's edge (fel-physics.md sec-source-filter). The angle is the knob and
! the grid-relative cut is what the kernel builds from, so the conversion happens here,
! where the beam that sets the angle and the grid that receives it both exist. Half the
! grid's Nyquist frequency is the angle lambda/(4 dx), which is what xcut = 1 means.

if (any(is_fel)) then
  block
    real(rp) th_mode, th_rho, th, th_p, margin, ang_per_xcut, ratio
    character(24) origin
    integer ie_first

    ! Genesis normalizes the shifted grid index by ngrid, so its x = 1 is the index ngrid
    ! itself, twice the Nyquist index, and xcut = 1 places the sigmoid's edge at the angle
    ! lambda/dx, off the grid. An edge at the angle theta is therefore xcut = theta dx /
    ! lambda. A first version of this line had lambda/(4 dx) and every edge landed four
    ! times wider than its name, which the Aramis benchmark hid and a second machine
    ! caught (FINDINGS 7.51). The grid comes from the input rather than from the built
    ! wavefront, since the conversion has to hold for every member of the field set.

    ie_first = findloc(is_fel, .true., dim = 1)
    ! A run whose field came from a file states no grid, so there is nothing to convert
    ! against. The filter is refused on that path below, and this stays finite.

    ang_per_xcut = 0
    if (run%winit%grid_half_width > 0 .and. run%winit%grid_n_pts > 1) &
        ang_per_xcut = fbeam%wavelength * (run%winit%grid_n_pts - 1) / (2 * run%winit%grid_half_width)

    call fel_filter_angles (fbeam, run%und_of(ie_first), fbeam%wavelength, th_mode, th_rho)

    ! The angle the stats split reports, whether or not the filter is on: it is the angle
    ! that separates the mode from the wide-angle emission of the point beamlets, and a
    ! run wants to see that separation most when it is not filtering.

    run%split_angle = max(th_mode, th_rho)
    run%split_origin = 'the mode angle'
    if (th_rho > th_mode) run%split_origin = 'the rho angle'
    if (run%global%source_filter_angle > 0) then
      run%split_angle = run%global%source_filter_angle
      run%split_origin = 'set by the deck'
    endif

    if (.not. run%global%source_filter) then
      continue
    else if (run%global%source_filter_xcut /= 0 .or. run%global%source_filter_ycut /= 0) then
      where (is_fel)
        run%und_of%filter%xcut = max(run%global%source_filter_xcut, 1.0e-12_rp)
        run%und_of%filter%ycut = max(run%global%source_filter_ycut, 1.0e-12_rp)
        run%und_of%filter%angle = 0
      end where
      call out_io (s_info$, r_name, 'Source filter: edge from xcut and ycut, ' // &
                   'validation-internal, at \es10.3\ rad on this grid.', &
                   'The cuts place the edge directly, so source_filter_tolerance is bypassed.', &
                   r_array = [run%global%source_filter_xcut * ang_per_xcut])

    else
      if (run%global%source_filter_angle > 0) then
        th = run%global%source_filter_angle
        origin = 'set by the deck'
      else if (max(th_mode, th_rho) <= 0) then
        call out_io (s_error$, r_name, 'SOURCE_FILTER IS ON AND ITS ANGLE CANNOT BE ' // &
                     'DERIVED: THE BEAM HAS NO TRANSVERSE SIZE OR NO CURRENT.', &
                     'POSSIBLE SOLUTION: SET SOURCE_FILTER_ANGLE.')
        err_flag = .true.;  return
      else

        ! The larger of the two angles is the one protected, since cutting into the mode
        ! loses real radiation while leaving artifact in only weakens the filter. Where
        ! the edge then goes is a separate question, and the tolerance answers it: the
        ! sigmoid's half-amplitude point on the protected angle passes a quarter of the
        ! intensity there, which is a filter taking three quarters of the source at the
        ! very angle it was told to keep. The edge is moved out until the transmission
        ! across that angle is what the tolerance allows.

        th_p = max(th_mode, th_rho)
        margin = fel_filter_margin (run%global%source_filter_tolerance, &
                                    run%global%source_filter_width)
        if (margin <= 0) then
          call out_io (s_error$, r_name, 'SOURCE_FILTER_WIDTH = \es12.3\ IS TOO LARGE FOR ' // &
                       'SOURCE_FILTER_TOLERANCE = \es12.3\ : NO EDGE HOLDS THAT TRANSMISSION.', &
                       'POSSIBLE SOLUTION: USE SOURCE_FILTER_WIDTH BELOW \es9.2\ , OR RAISE THE TOLERANCE.', &
                       r_array = [run%global%source_filter_width, &
                                  run%global%source_filter_tolerance, &
                                  run%global%source_filter_width / (1 - margin)])
          err_flag = .true.;  return
        endif
        th = th_p / margin
        origin = run%split_origin
      endif

      ! Their ratio goes as sqrt(z_R/L_g), 15 on the Aramis benchmark and near 1 on a
      ! diffraction-dominated machine, and a case far outside that range is unlike the one
      ! the default was measured on (doc/startup-noise.md).

      ratio = 0
      if (th_rho > 0) ratio = th_mode / th_rho
      call out_io (s_info$, r_name, 'Source filter: mode angle \es10.3\ rad, ' // &
                   'rho angle \es10.3\ rad, ratio \f8.2\ .', &
                   r_array = [th_mode, th_rho, ratio])
      if (run%global%source_filter_angle > 0) then
        call out_io (s_info$, r_name, 'Edge at \es12.5\ rad, set by the deck, so ' // &
                     'source_filter_tolerance is bypassed. Sigmoid on axis \f7.4\ .', &
                     r_array = [th, 1 / (1 + exp(-1 / run%global%source_filter_width))])
      else
        call out_io (s_info$, r_name, 'Protected angle \es12.5\ rad, ' // trim(origin) // &
                     '. Tolerance \f7.4\ puts the edge at \es12.5\ rad. Sigmoid on axis \f7.4\ .', &
                     r_array = [th_p, run%global%source_filter_tolerance, th, &
                                1 / (1 + exp(-1 / run%global%source_filter_width))])
      endif
      if (ratio > 0 .and. (ratio < 0.3_rp .or. ratio > 3.0_rp)) then
        call out_io (s_warn$, r_name, 'The two filter angles differ by more than the range ' // &
                     'the default was measured over (0.3 to 3).', &
                     'Check the power inside the edge against the total before trusting the run.')
      endif

      where (is_fel)
        run%und_of%filter%xcut = th / ang_per_xcut
        run%und_of%filter%ycut = th / ang_per_xcut
        run%und_of%filter%angle = th
      end where
    endif
  end block
endif

! More than one slice means a time-dependent run with slippage active. One slice is the
! steady state and fel_apply_slippage is a no-op.

timerun = (nslice > 1)
do ih = 1, n_harm       ! Same values for every field: one window, lockstep rotation
                        ! in fundamental-wavelength units (Genesis's one Control::sample).
  ffield(ih)%slip%timerun = timerun
  ffield(ih)%slip%sample = fbeam%n_wavelength
enddo
run%gamma0_ref = fel_gamma0(fbeam)
gamma0_ref = run%gamma0_ref

call check_wake_window ()
if (err_flag) return
run%any_unavg = any(fel_mode == unaveraged$ .and. is_fel)

! The collective terms are not wired into the unaveraged step (fel-physics.md
! sec-unaveraged), and a mixed line would apply them in some segments and silently
! drop them in others. Refuse.

if (run%any_unavg .and. (wake_on .or. sc_nz >= 1 .or. sc_longrange)) then
  call out_io (s_error$, r_name, 'WAKES/SPACE CHARGE ARE NOT WIRED INTO THE UNAVERAGED MODE', &
                                 '(SEE fel-physics.md sec-unaveraged).', &
                                 'POSSIBLE SOLUTION: TURN THEM OFF.')
  err_flag = .true.;  return
endif


! The rest of the schedule: drift autophasing. Interludes accumulate Lz, and the last
! interlude before each undulator gets floor(Lz/(2*gamma0^2*lambda)) + 1 wavelengths
! (Lattice.cpp:171-174, guarded there by Lz > 0). The end-of-lattice fixup
! (Lattice.cpp:191-193) is unguarded in Genesis: the last element always gets
! floor(Lz/(2*gamma0^2*lambda)) + 1, which is +1 even with no trailing interlude at all
! ("autophasing is applied in case for [a] second, succeeding run"). Transcribed as is:
! omitting that +1 leaves the field record one rotation short at the very end, found the
! hard way against the single-segment time-dependent run. (Citations kept at the lines:
! this quirk's exactness matters here, at the call site. Manual sec-slippage.)

allocate (run%ele_slip(branch%n_ele_track))
allocate (run%fel_zoff(branch%n_ele_track), run%light_corr(branch%n_ele_track))
ele_slip => run%ele_slip;  fel_zoff => run%fel_zoff;  light_corr => run%light_corr
ele_slip = 0
fel_zoff = 0
light_corr = 0
phase_rate = twopi / (2 * gamma0_ref**2 * wf%wavelength)
run%phase_rate = phase_rate
Lz = 0
prev_ie = 0

do ie = 1, branch%n_ele_track
  ele => branch%ele(ie)
  if (ele%value(l$) == 0) cycle
  if (is_fel(ie)) then
    if (Lz > 0 .and. prev_ie > 0) then
      ele_slip(prev_ie) = ele_slip(prev_ie) + floor(Lz / (2 * gamma0_ref**2 * wf%wavelength)) + 1
    endif

    ! The off-phase knob (fel-physics.md sec-phasing): the wiggler's own z_offset, standard
    ! Bmad misalignment (girder-composed _tot form). Anchored at the nominal position:
    ! the entry phase shifts by the displaced upstream break, the exit unshifts for the
    ! displaced downstream one, so everything downstream stays anchored. The knob must
    ! fit inside its breaks -- a z_offset that walks the element out of its gap is
    ! geometry, not phasing.

    fel_zoff(ie) = ele%value(z_offset_tot$)
    if (fel_zoff(ie) /= 0 .and. prev_ie > 0 .and. abs(fel_zoff(ie)) >= Lz .and. Lz > 0) then
      call out_io (s_error$, r_name, 'THE Z_OFFSET OF FEL ELEMENT ' // trim(ele%name), &
                   '(\es10.2\ m) EXCEEDS ITS UPSTREAM BREAK (\es10.2\ m).', r_array = [fel_zoff(ie), Lz])
      err_flag = .true.;  return
    endif
    if (fel_zoff(ie) /= 0 .and. prev_ie == 0) then
      call out_io (s_error$, r_name, 'THE Z_OFFSET KNOB ON THE FIRST ELEMENT ' // trim(ele%name), &
                   'HAS NO UPSTREAM BREAK TO DISPLACE INTO; GIVE THE LATTICE A LEADING BREAK.')
      err_flag = .true.;  return
    endif
    Lz = 0
  else
    Lz = Lz + ele%value(l$)
  endif
  prev_ie = ie
enddo
if (prev_ie > 0) then
  ele_slip(prev_ie) = ele_slip(prev_ie) + floor(Lz / (2 * gamma0_ref**2 * wf%wavelength)) + 1
endif

call setup_break_geometry ()   ! Chicane breaks: chord vs arc from ele%floor, the
                               ! delay's rotations, the light-path correction, and
                               ! the closed-bump and genesis-model refusals.
if (err_flag) return

! Diagnostics file, one row per slice per record at Genesis's record positions, slices in
! time-window order.

! The tracking window (global%track_start/track_end: Tao's names, with Genesis zstop
! parity). Resolved through Bmad's own locator. The schedule above was built on the
! full lattice, so a windowed run composes exactly with the full one. The walk
! simply covers [i_start, i_end], and no end-of-lattice fixup moves.

run%i_start = 1
run%i_end = branch%n_ele_track
if (run%global%track_start /= '') then
  call resolve_window_ele (run%global%track_start, 'track_start', run%i_start)
  if (err_flag) return
endif
if (run%global%track_end /= '') then
  call resolve_window_ele (run%global%track_end, 'track_end', run%i_end)
  if (err_flag) return
endif
if (run%i_start > run%i_end) then
  call out_io (s_error$, r_name, 'TRACK_START (ELEMENT \i0\ ) IS PAST TRACK_END (ELEMENT \i0\ ).', &
               i_array = [run%i_start, run%i_end])
  err_flag = .true.;  return
endif

! The slices a frame carries. The default pair is the whole window, and -1 as the last
! reads as the last slice so a range needs one convention and not two. The dumps that a
! run restarts from are whole whatever this says, so the range is refused only against
! the window it is a range of.

run%dump_is1 = run%global%dump_slice_first
run%dump_is2 = run%global%dump_slice_last
if (run%dump_is2 == -1) run%dump_is2 = run%nslice
if (run%dump_is1 < 1 .or. run%dump_is2 > run%nslice) then
  call out_io (s_error$, r_name, 'DUMP_SLICE_FIRST \i0\ AND DUMP_SLICE_LAST \i0\ ARE NOT ' // &
               'INSIDE THE WINDOW OF \i0\ SLICES.', &
               i_array = [run%dump_is1, run%dump_is2, run%nslice])
  err_flag = .true.;  return
endif
if (run%dump_is1 > run%dump_is2) then
  call out_io (s_error$, r_name, 'DUMP_SLICE_FIRST \i0\ IS PAST DUMP_SLICE_LAST \i0\ .', &
               i_array = [run%dump_is1, run%dump_is2])
  err_flag = .true.;  return
endif

call setup_diagnostics ()
if (err_flag) return

! The FP32 lockstep instrument. It carries a twin for each advance, so a lattice of
! either mode or of both is instrumented, each segment by the twin that mirrors it. The
! unaveraged twin prices that mode's particle path, where its reformulations live, and a
! field row of its own. doc/validation.md records both.
!
! Freerun compounds a single-precision state across steps and the unaveraged twin does
! not carry one, so that combination is refused rather than reported as a lockstep.

if (run%global%fp32_check == 'freerun' .and. run%any_unavg) then
  call out_io (s_error$, r_name, 'FP32_CHECK = "freerun" DOES NOT COVER THE UNAVERAGED MODE.', &
               'THE UNAVERAGED TWIN IS A LOCKSTEP: IT IS REBUILT FROM THE FP64 STATE EVERY', &
               'RECORD STEP AND CARRIES NOTHING ACROSS THEM. USE FP32_CHECK = "lockstep".')
  err_flag = .true.;  return
endif
call fel_fp32_setup (run%fp32, run%global%fp32_check, run%global%fp32_mutate, run%nslice, &
                     fel_p0_mc(fbeam) / fel_gamma0(fbeam) * fbeam%slice_spacing, trim(out_root), err_flag)
if (err_flag) return

! The device backend. Everything its kernels do not cover is refused here --
! an unsupported configuration stops the run and never quietly takes the CPU path.
! Element wakes are refused at the element (only the walk sees them), and the grid
! refusal, naming the nearest supported size, comes back from the backend itself.
! The field set is covered whole: every harmonic member and both polarization planes
! ride the device, and their combination is refused above for the CPU and the device
! alike. Slice migration is covered without a kernel: it runs on the host at an
! element's last step after the walk reads the beam back (do_migrate), and the next
! element re-uploads into a rectangle grown to the new largest fill.

if (run%global%device /= '' .and. run%global%device /= 'off') then
  if (run%any_unavg) then
    call out_io (s_error$, r_name, 'DEVICE = "' // trim(run%global%device) // '" DOES NOT COVER THE UNAVERAGED MODE.')
    err_flag = .true.;  return
  endif
  if (run%global%source_model == 'coherent') then
    call out_io (s_error$, r_name, 'DEVICE = "' // trim(run%global%device) // '" DOES NOT COVER THE COHERENT SOURCE MODEL.')
    err_flag = .true.;  return
  endif
  if (run%chamber_wake%on) then
    call out_io (s_error$, r_name, 'DEVICE = "' // trim(run%global%device) // '" DOES NOT COVER WAKES.')
    err_flag = .true.;  return
  endif
  if (run%global%source_filter_mutate) then
    call out_io (s_error$, r_name, 'DEVICE = "' // trim(run%global%device) // &
                 '" DOES NOT COVER SOURCE_FILTER_MUTATE, WHICH IS THE CPU CHECK''S OWN HOOK.')
    err_flag = .true.;  return
  endif
  if (any(run%sc_here)) then
    call out_io (s_error$, r_name, 'DEVICE = "' // trim(run%global%device) // '" DOES NOT COVER SPACE CHARGE.')
    err_flag = .true.;  return
  endif
  if (bmad_com%radiation_damping_on .or. bmad_com%radiation_fluctuations_on) then
    call out_io (s_error$, r_name, 'DEVICE = "' // trim(run%global%device) // '" DOES NOT COVER SPONTANEOUS RADIATION.')
    err_flag = .true.;  return
  endif
  if (run%global%keep_escaped_field) then
    call out_io (s_error$, r_name, 'DEVICE = "' // trim(run%global%device) // '" DOES NOT COVER THE ESCAPED-FIELD BANK.')
    err_flag = .true.;  return
  endif
  if (run%global%transport_model /= 'bmad') then
    call out_io (s_error$, r_name, 'DEVICE = "' // trim(run%global%device) // '" CARRIES THE BMAD TRANSVERSE MAPS ONLY;', &
                                   'TRANSPORT_MODEL = "' // trim(run%global%transport_model) // '" IS VALIDATION-INTERNAL.')
    err_flag = .true.;  return
  endif
  block
    integer ngrid_dev(3)
    ngrid_dev = wavefront_shape(run%ffield(1)%wf)
    call fel_device_setup (run%dev, run%global%device, fbeam, ngrid_dev(1), &
                           run%ffield(1)%slip%sample, run%ffield(1:run%n_harm)%harm, &
                           merge(2, 1, run%two_pol), run%fp32%iu, &
                           run%global%device_timing, run%global%device_dep_mutate, err_flag)
  end block
  if (err_flag) return
endif

!------------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!+
! Subroutine resolve_window_ele (locator, which, ix)
!
! Routine to resolve one tracking-window locator to a single tracked-element index,
! refused when it matches nothing or more than one element.
!-

subroutine resolve_window_ele (locator, which, ix)

type (ele_pointer_struct), allocatable :: eles(:)
character(*) locator, which
integer ix, n_loc
logical lerr

call lat_ele_locator (locator, lat, eles, n_loc, lerr)
if (lerr .or. n_loc == 0) then
  call out_io (s_error$, r_name, upcase(trim(which)) // ' MATCHES NO ELEMENT: ' // trim(locator))
  err_flag = .true.;  return
endif
if (n_loc > 1) then
  call out_io (s_error$, r_name, upcase(trim(which)) // ' MATCHES MORE THAN ONE ELEMENT: ' // trim(locator))
  err_flag = .true.;  return
endif
if (eles(1)%ele%ix_ele < 1 .or. eles(1)%ele%ix_ele > branch%n_ele_track) then
  call out_io (s_error$, r_name, upcase(trim(which)) // ' IS NOT A TRACKED ELEMENT: ' // trim(locator))
  err_flag = .true.;  return
endif
ix = eles(1)%ele%ix_ele

end subroutine resolve_window_ele

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine check_wake_window ()
!
! Routine to check the Bmad element wakes against the time window. Element sr wakes act
! across the whole window (fel-physics.md sec-seamwake): all slices concatenate into one bunch
! in global window coordinates and Bmad's wake machinery applies unmodified. What is
! checked here: lr (multi-bunch) wakes are not supported. A pseudomode wake
! whose z_max is shorter than the window would have Bmad kill the bunch mid-run. A
! z_long table narrower than the window would overflow its binning grid the same way.
! Runs after the beam is built (the window length is the subject). setup_fel_elements
! runs before it (two_pol must precede the field).
!-

subroutine check_wake_window ()

type (ele_struct), pointer :: w
integer je

do je = 1, branch%n_ele_track
  w => pointer_to_wake_ele(branch%ele(je))
  if (.not. associated(w)) cycle
  if (allocated(w%wake%lr%mode)) then
    if (size(w%wake%lr%mode) > 0) then
      call out_io (s_error$, r_name, 'LR (MULTI-BUNCH) WAKES ARE NOT SUPPORTED;', &
                   'REMOVE THEM FROM: ' // trim(w%name))
      err_flag = .true.;  return
    endif
  endif
  if (w%wake%sr%z_max > 0 .and. size(fbeam%slice) * fbeam%slice_spacing > w%wake%sr%z_max) then
    call out_io (s_error$, r_name, 'THE TIME WINDOW IS LONGER THAN THIS ELEMENT''S SR WAKE', &
                 'Z_MAX CAN HANDLE: ' // trim(w%name))
    err_flag = .true.;  return
  endif
  if (w%wake%sr%z_long%dz > 0 .and. &
      size(fbeam%slice) * fbeam%slice_spacing > w%wake%sr%z_long%z0) then
    call out_io (s_error$, r_name, 'THE TIME WINDOW IS LONGER THAN THIS ELEMENT''S Z_LONG', &
                 'WAKE TABLE EXTENT Z0: ' // trim(w%name))
    err_flag = .true.;  return
  endif
enddo

end subroutine check_wake_window

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine setup_break_geometry ()
!
! Routine to set up chicane breaks (fel-physics.md sec-phasing): a break whose elements bend the
! reference (sbends, patches) detours the beam while the radiation goes straight. The
! light's path is the chord between the flanking undulator faces, from ele%floor,
! never the reference arc that vec(5) is measured against. The arc-minus-chord
! delay is charged as whole-wavelength window rotations (Genesis's chicane
! semantics: "always autophasing") on the break's last element, which also takes
! the light-path drift correction. Absolute mode adds the delay's carrier phase in
! the walk. Only a closed bump keeps the light on the next undulator's axis.
! Anything else is refused, as is any geometry element under the
! genesis-model interludes (Genesis's drift/quad set cannot represent it).
!-

subroutine setup_break_geometry ()

real(rp) arc
integer i0, ie_g, last_in_break
logical geom

!

i0 = 0                    ! Break start: exit face of the last FEL element (0 = origin).
arc = 0;  geom = .false.;  last_in_break = 0

do ie_g = 1, branch%n_ele_track
  ele => branch%ele(ie_g)
  if (is_fel(ie_g)) then
    if (geom .and. last_in_break > 0) then
      call close_geometry_break (i0, last_in_break, arc)
      if (err_flag) return
    endif
    i0 = ie_g;  arc = 0;  geom = .false.;  last_in_break = 0
  elseif (ele%value(l$) /= 0 .or. ele%key == patch$) then
    arc = arc + ele%value(l$)
    last_in_break = ie_g
    if (ele%key == sbend$ .or. ele%key == patch$) then
      geom = .true.
      if (interlude_model == 'genesis') then
        call out_io (s_error$, r_name, 'GEOMETRY ELEMENT ' // trim(ele%name) // ' (A BEND OR PATCH) INSIDE A', &
                     'GENESIS-MODEL INTERLUDE: GENESIS''S DRIFT/QUAD SET CANNOT REPRESENT IT.', &
                     'POSSIBLE SOLUTION: USE INTERLUDE_MODEL = "bmad" (THE SEAM TRACKS IT EXACTLY).')
        err_flag = .true.;  return
      endif
    endif
  endif
enddo
! A trailing geometry break (no following undulator) needs no phasing: there is no
! next segment to phase against.

end subroutine setup_break_geometry

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine close_geometry_break (ib0, iblast, arc)
!
! Routine to close one geometry break: the chord, the closed-bump refusal, the delay's
! rotations and the light-path correction (see setup_break_geometry's header).
!-

subroutine close_geometry_break (ib0, iblast, arc)

type (ele_struct), pointer :: e1, e2
real(rp) arc, chord, dvec(3), axis(3), delay_geo, tol
integer ib0, iblast

!

e1 => branch%ele(ib0)         ! Exit face of the upstream FEL element (or the origin).
e2 => branch%ele(iblast)      ! Exit face of the break's last element = the entry face
                              ! of the next FEL element.
dvec = e2%floor%r - e1%floor%r
chord = norm2(dvec)
axis = e1%floor%w(:,3)
tol = 1e-9_rp * max(1.0_rp, chord)

if (abs(e2%floor%theta - e1%floor%theta) > 1e-9_rp .or. abs(e2%floor%phi - e1%floor%phi) > 1e-9_rp .or. &
    abs(e2%floor%psi - e1%floor%psi) > 1e-9_rp .or. norm2(dvec - chord * axis) > tol) then
  call out_io (s_error$, r_name, 'THE BREAK ENDING AT ' // trim(branch%ele(iblast)%name) // ' BENDS THE REFERENCE', &
               'AND IS NOT A CLOSED BUMP: THE RADIATION WOULD LEAVE THE NEXT UNDULATOR''S AXIS.', &
               'ONLY CLOSED-BUMP CHICANES ARE MODELED; CLOSE THE GEOMETRY OR STRAIGHTEN THE LINE.')
  err_flag = .true.;  return
endif

delay_geo = arc - chord
ele_slip(iblast) = ele_slip(iblast) + floor(delay_geo / wf%wavelength)
light_corr(iblast) = delay_geo

end subroutine close_geometry_break

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine setup_diagnostics ()
!
! Routine to set up the diagnostics (fel-physics.md sec-stats): resolve the dump-at lists through
! Bmad's own lat_ele_locator (class::name syntax for free) and precompute the exact
! record and element-end counts by replaying the walk's skip rule, so the stats arrays
! are sized once, never grown. An entry matching nothing is refused.
!-

subroutine setup_diagnostics ()

type (ele_pointer_struct), allocatable :: eles(:)
integer i, j, n_loc, nstep_r, istep_r
real(rp) comb_r, z_r, z_last_r, dz_r
logical derr

!

allocate (run%dump_beam_here(0:branch%n_ele_track), run%dump_field_here(0:branch%n_ele_track))
dump_beam_here => run%dump_beam_here;  dump_field_here => run%dump_field_here
dump_beam_here = .false.;  dump_field_here = .false.

do i = 1, size(dump_beam_at)
  if (dump_beam_at(i) == '') cycle
  call lat_ele_locator (dump_beam_at(i), lat, eles, n_loc, derr)
  if (derr .or. n_loc == 0) then
    call out_io (s_error$, r_name, 'DUMP_BEAM_AT ENTRY MATCHES NO ELEMENT: ' // trim(dump_beam_at(i)))
    err_flag = .true.;  return
  endif
  do j = 1, n_loc
    if (eles(j)%ele%ix_ele >= 1 .and. eles(j)%ele%ix_ele <= branch%n_ele_track) &
                            dump_beam_here(eles(j)%ele%ix_ele) = .true.
  enddo
enddo

do i = 1, size(dump_field_at)
  if (dump_field_at(i) == '') cycle
  call lat_ele_locator (dump_field_at(i), lat, eles, n_loc, derr)
  if (derr .or. n_loc == 0) then
    call out_io (s_error$, r_name, 'DUMP_FIELD_AT ENTRY MATCHES NO ELEMENT: ' // trim(dump_field_at(i)))
    err_flag = .true.;  return
  endif
  do j = 1, n_loc
    if (eles(j)%ele%ix_ele >= 1 .and. eles(j)%ele%ix_ele <= branch%n_ele_track) &
                            dump_field_here(eles(j)%ele%ix_ele) = .true.
  enddo
enddo

! The record count replays the walk's skip rule, its window, its z arithmetic and
! the comb rule (fel_comb_take, the one authority), so the stats arrays are
! exact-sized in every mode -- never grown, never padded.

comb_r = run%global%comb_ds_save
z_r = branch%ele(run%i_start - 1)%s
z_last_r = -1e30_rp
nrec_stats = 0
nend_stats = 0
if (fel_comb_take(comb_r, z_r, z_last_r, .false.)) nrec_stats = nrec_stats + 1
do i = run%i_start, run%i_end
  ele => branch%ele(i)
  wake_src => pointer_to_wake_ele(ele)
  if (ele%value(l$) == 0 .and. .not. associated(wake_src)) cycle
  nend_stats = nend_stats + 1
  if (is_fel(i)) then
    nstep_r = max(1, nint(ele%value(num_steps$)))
    dz_r = ele%value(l$) / nstep_r
    do istep_r = 1, nstep_r
      z_r = z_r + dz_r
      if (fel_comb_take(comb_r, z_r, z_last_r, istep_r == nstep_r)) nrec_stats = nrec_stats + 1
    enddo
  else
    z_r = z_r + ele%value(l$)
    if (fel_comb_take(comb_r, z_r, z_last_r, .true.)) nrec_stats = nrec_stats + 1
  endif
enddo

run%nrec_stats = nrec_stats;  run%nend_stats = nend_stats
call fel_stats_init (stats, nrec_stats, nend_stats, nslice, fbeam%p0c, two_pol, harmonics(2:n_harm))
stats%split_angle = run%split_angle

! The field's reduced projections, where the run asked for them. One entry per member and
! plane, named as the recorder fills them.

if (run%global%dump_reduced) then
  block
    character(12), allocatable :: rname(:)
    integer nred, k, ih, ngrid(3)
    nred = n_harm * merge(2, 1, two_pol)
    allocate (rname(nred))
    k = 0
    do ih = 1, n_harm
      k = k + 1
      write (rname(k), '(a, i0, a)') 'h', harmonics(ih), '_x'
      if (two_pol) then
        k = k + 1
        write (rname(k), '(a, i0, a)') 'h', harmonics(ih), '_y'
      endif
    enddo
    ! The grid is the built wavefront's and not the deck's: a derived grid or one read
    ! from a field file is what the projections are summed over.
    ngrid = wavefront_shape(run%ffield(1)%wf)
    call fel_stats_reduced_init (stats, ngrid(1), ngrid(2), rname)
  end block
endif
allocate (run%bdiag_arr(nslice), run%fpow_arr(nslice), run%fonax_arr(nslice))

end subroutine setup_diagnostics

end subroutine fel_setup_schedule


end module fel_setup_mod
