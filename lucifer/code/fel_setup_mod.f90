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
logical, save, private :: fel_attributes_registered = .false.

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
! FEL elements carry tracking_method = custom, and Bmad's bookkeeping (the reference
! time/energy pass inside bmad_parser, any track1 at the seam) resolves custom tracking
! through track1_custom_ptr. Point it at the standard periodic-wiggler kernel, so the
! element behaves as the plain Bmad wiggler it is everywhere EXCEPT inside this
! driver's own FEL walk. In particular the reference time acquires the resonant
! undulation delay from Bmad's own code, not from anything written here.

! Both hooks are needed: mat6_calc_method resolves to custom too (auto follows the
! tracking method), and make_mat6 calls through a null make_mat6_custom_ptr otherwise.

track1_custom_ptr => fel_ele_as_wiggler
make_mat6_custom_ptr => fel_mat6_as_wiggler

! The FEL mode and unaveraged parameters live on the LATTICE (manual sec:element),
! registered program-side so no lattice declares them. The same slot index serves
! wigglers and undulators.

if (.not. fel_attributes_registered) then
  call set_custom_attribute_name ('WIGGLER::FEL_TRACKING', err, 1)
  if (.not. err) call set_custom_attribute_name ('UNDULATOR::FEL_TRACKING', err, 1)
  if (.not. err) call set_custom_attribute_name ('WIGGLER::FEL_STEPS_PER_PERIOD', err, 2)
  if (.not. err) call set_custom_attribute_name ('UNDULATOR::FEL_STEPS_PER_PERIOD', err, 2)
  if (.not. err) call set_custom_attribute_name ('WIGGLER::FEL_RAMP_PERIODS', err, 3)
  if (.not. err) call set_custom_attribute_name ('UNDULATOR::FEL_RAMP_PERIODS', err, 3)
  if (err) then
    call out_io (s_error$, r_name, 'COULD NOT REGISTER THE FEL LATTICE ATTRIBUTES.')
    err_flag = .true.;  return
  endif
  fel_attributes_registered = .true.
endif

! err_flag matters: bmad_parser reports attribute errors (e.g. a wake on an element
! type that cannot carry one) and RETURNS. Without the check the run continues on a
! partial lattice, found when an example's drift wakes silently never attached.

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

! The field set (manual sec:field-set). Validate the harmonics request, then allocate
! the set and point the fundamental aliases at entry 1 BEFORE any field construction:
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
if (n_harm > 1 .and. any(fel_mode == fel_unaveraged$ .and. is_fel)) then
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

! The source model (manual sec:coherent-source): validated, then stamped onto every
! FEL element. v1 scope refusals, each by name: the coherent source carries the
! FUNDAMENTAL of one polarization, and it cannot live inside the unaveraged referee
! (variance reduction inside the explicit-everything mode). Even harmonics are invalid
! in the method itself -- F(z,0) = 0 -- and odd ones are a named follow-on.

select case (run%global%source_model)
case ('deposit')
case ('coherent')
  if (any(fel_mode == fel_unaveraged$ .and. is_fel)) then
    call out_io (s_error$, r_name, 'SOURCE_MODEL = "coherent" WITH AN UNAVERAGED ELEMENT: THE', &
                                   'UNAVERAGED MODE IS THE EXPLICIT REFEREE; VARIANCE REDUCTION INSIDE IT IS REFUSED.')
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

    ! MEASURED, not assumed (check_coherent's SASE experiment): the coherent source
    ! understates dark-start startup by ~175x on the reference configuration.
    ! Spontaneous, spatially-INCOHERENT emission dominates SASE startup and is
    ! exactly what the coherent model drops. The slice bunch factor's physical noise
    ! (Fawley, <|B|^2> N_lambda = 1) is present but is not the dominant seed.
    ! Refused by name. Seeded runs are fully supported.

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
! segment is a wiggler/undulator element with tracking_method = custom: Bmad's own
! semantics for program-supplied tracking, which this driver is. The wiggler sanity
! assertions are enforced: a wiggler with zero b_max or l_period would silently get
! factor = 0 in Bmad's own kernel (no resonance, no error), and a fieldmap field_calc
! gets osc_amplitude without focusing. Both are refused by name. The stored k1x/k1y
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

do je = 1, branch%n_ele_track
  w => branch%ele(je)
  if (.not. (w%key == wiggler$ .or. w%key == undulator$)) cycle
  if (w%tracking_method /= custom$) cycle

  ! The FEL mode and unaveraged parameters, from the element's own attributes.
  ! fel_tracking: unset/0 = averaged with the bmad_standard kernel's transverse maps
  ! (Bmad's own kernel is the default). 1 = unaveraged. -1 = averaged with the
  ! transcribed-Genesis maps (validation-internal: the Genesis tiers require
  ! transcription-level transport, and no production lattice writes it).

  rv = value_of_attribute(w, 'FEL_TRACKING', err_a)
  if (err_a) then
    err_flag = .true.;  return
  endif
  fel_mode(je) = nint(rv)
  if (abs(rv - fel_mode(je)) > 1e-9_rp .or. fel_mode(je) < fel_transcribed$ .or. &
      fel_mode(je) > fel_unaveraged$) then
    call out_io (s_error$, r_name, 'FEL_TRACKING MUST BE -1 (TRANSCRIBED MAPS, VALIDATION-INTERNAL),', &
                 '0/UNSET (AVERAGED, BMAD_STANDARD KERNEL MAPS) OR 1 (UNAVERAGED), AT ELEMENT: ' // trim(w%name))
    err_flag = .true.;  return
  endif

  rv = value_of_attribute(w, 'FEL_STEPS_PER_PERIOD', err_a)
  if (err_a) then
    err_flag = .true.;  return
  endif
  fel_spp(je) = nint(rv)
  if (fel_spp(je) == 0) fel_spp(je) = 20
  if (fel_spp(je) < 10) then
    call out_io (s_error$, r_name, 'FEL_STEPS_PER_PERIOD IS BELOW THE FLOOR OF 10 (MINERVA''S ENVELOPE),', &
                 'AT ELEMENT: ' // trim(w%name))
    err_flag = .true.;  return
  endif

  ! fel_ramp_periods: an attribute's unset value is 0, and a silent hard edge would
  ! reintroduce the K/gamma handoff hazard by omission. Thus unset/0 means the default
  ! of 2 periods, and a TRUE hard edge (the mutation/test configuration) must be asked
  ! for by name with the explicit sentinel -1.

  rv = value_of_attribute(w, 'FEL_RAMP_PERIODS', err_a)
  if (err_a) then
    err_flag = .true.;  return
  endif
  if (rv == 0) then
    fel_ramp(je) = 2
  elseif (rv == -1) then
    fel_ramp(je) = 0
  elseif (rv < 0) then
    call out_io (s_error$, r_name, 'FEL_RAMP_PERIODS MUST BE POSITIVE, 0/UNSET (DEFAULT 2), OR THE', &
                 'HARD-EDGE TEST SENTINEL -1, AT ELEMENT: ' // trim(w%name))
    err_flag = .true.;  return
  else
    fel_ramp(je) = rv
  endif

  ! The wiggler sanity assertions live in fel_assert_wiggler_sane, ONE authority. It is
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
  ! helicity, scaled by ku^2 as Genesis's unroll does (manual sec:element).

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
    if (fel_mode(je) == fel_transcribed$) then
      call out_io (s_error$, r_name, 'THE TRANSCRIBED-GENESIS MAPS (FEL_TRACKING = -1) KNOW NO TILT', &
                   '(GENESIS HAS NONE); USE THE DEFAULT MAPS ON: ' // trim(w%name))
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
  call out_io (s_error$, r_name, 'THE LATTICE HAS NO FEL ELEMENTS (WIGGLER/UNDULATOR WITH TRACKING_METHOD = CUSTOM).', &
               'POSSIBLE SOLUTION: SET REFERENCE_RUN = T FOR A DELIBERATE NO-FEL RUN (BMAD TRACKS EVERYTHING).')
  err_flag = .true.;  return
endif

end subroutine setup_fel_elements

end subroutine fel_setup_lattice

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_setup_schedule (run, err_flag)
!
! Routine to build everything that needs BOTH the lattice and the built starting state:
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
integer sc_ngrid, sc_nz, sc_nphi
real(rp) gamma0_ref, phase_rate, Lz
integer nslice, n_harm, ie, ih, prev_ie, nrec_stats, nend_stats
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
wake_on = run%wake_init%on
wake_loss = run%wake_init%loss
wake_radius = run%wake_init%radius
wake_conductivity = run%wake_init%conductivity
wake_relaxation = run%wake_init%relaxation
wake_roundpipe = run%wake_init%roundpipe
wake_material = run%wake_init%material
wake_gap = run%wake_init%gap
wake_lgap = run%wake_init%lgap
wake_hrough = run%wake_init%hrough
wake_lrough = run%wake_init%lrough
sc_rmax = run%sc_init%rmax
sc_ngrid = run%sc_init%ngrid
sc_nz = run%sc_init%nz
sc_nphi = run%sc_init%nphi
sc_longrange = run%sc_init%longrange


coll%efield%on = (sc_nz >= 1 .or. sc_longrange)
coll%efield%rmax = sc_rmax
coll%efield%ngrid = sc_ngrid
coll%efield%nz = sc_nz
coll%efield%nphi = sc_nphi
coll%efield%longrange = sc_longrange

if (coll%efield%on .and. interlude_model == 'bmad') then
  call out_io (s_info$, r_name, 'Note: space charge acts inside undulators and genesis-model', &
               'interludes only; the Bmad seam''s own space charge is a named follow-on.')
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
  ! convolution. NOTE the s = 0 entries carry the Bane self-slice half factor
  ! (fel_wake_init halves them). A plain W(z) table wants the unhalved value.

  if (run%wake_init%write_kernels /= '') then
    block
      integer iu_k, i_k
      open (newunit = iu_k, file = trim(run%wake_init%write_kernels), action = 'write')
      write (iu_k, '(a)') '# s [m]   wakeres   wakegeo   wakerou   [eV/(m electron)]; s=0 rows are HALVED (Bane self-slice)'
      do i_k = 1, coll%wake%ns
        write (iu_k, '(4es24.15e3)') coll%wake%ds * (i_k-1), coll%wake%wakeres(i_k), &
                                     coll%wake%wakegeo(i_k), coll%wake%wakerou(i_k)
      enddo
      close (iu_k)
      call out_io (s_info$, r_name, 'Wrote wake kernels: ' // trim(run%wake_init%write_kernels))
    end block
  endif
endif

! More than one slice means a time-dependent run with slippage active. One slice is the
! steady state and fel_apply_slippage is a no-op.

timerun = (nslice > 1)
do ih = 1, n_harm       ! Same values for every field: one window, lockstep rotation
                        ! in fundamental-wavelength units (Genesis's one Control::sample).
  ffield(ih)%slip%timerun = timerun
  ffield(ih)%slip%sample = fbeam%slice_spacing / fbeam%wavelength
enddo
run%gamma0_ref = fel_gamma0(fbeam)
gamma0_ref = run%gamma0_ref

call check_wake_window ()
if (err_flag) return
run%any_unavg = any(fel_mode == fel_unaveraged$ .and. is_fel)

! The collective terms are not wired into the unaveraged step (fel-physics.md
! sec:unaveraged), and a mixed line would apply them in some segments and silently
! drop them in others. Refuse by name.

if (run%any_unavg .and. (wake_on .or. sc_nz >= 1 .or. sc_longrange)) then
  call out_io (s_error$, r_name, 'WAKES/SPACE CHARGE ARE NOT WIRED INTO THE UNAVERAGED MODE', &
                                 '(SEE fel-physics.md sec-unaveraged).', &
                                 'POSSIBLE SOLUTION: TURN THEM OFF.')
  err_flag = .true.;  return
endif


! The rest of the schedule: drift autophasing. Interludes accumulate Lz, and the last
! interlude before each undulator gets floor(Lz/(2*gamma0^2*lambda)) + 1 wavelengths
! (Lattice.cpp:171-174, guarded there by Lz > 0). The end-of-lattice fixup
! (Lattice.cpp:191-193) is UNGUARDED in Genesis: the last element always gets
! floor(Lz/(2*gamma0^2*lambda)) + 1, which is +1 even with no trailing interlude at all
! ("autophasing is applied in case for [a] second, succeeding run"). Transcribed as is:
! omitting that +1 leaves the field record one rotation short at the very end, found the
! hard way against the single-segment time-dependent run. (Citations kept at the lines:
! this quirk's exactness matters here, at the call site. Manual sec:slippage.)

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

    ! The off-phase knob (manual sec:phasing): the wiggler's own z_offset, standard
    ! Bmad misalignment (girder-composed _tot form). Anchored at the NOMINAL position:
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
! parity). Resolved through Bmad's own locator. The SCHEDULE above was built on the
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

call setup_diagnostics ()
if (err_flag) return

!------------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
!+
! Subroutine resolve_window_ele (locator, which, ix)
!
! Routine to resolve one tracking-window locator to a single tracked-element index,
! refused by name when it matches nothing or more than one element.
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
! across the WHOLE window (manual sec:seamwake): all slices concatenate into one bunch
! in global window coordinates and Bmad's wake machinery applies unmodified. What is
! checked here, by name: lr (multi-bunch) wakes are not supported. A pseudomode wake
! whose z_max is shorter than the window would have Bmad kill the bunch mid-run. A
! z_long table narrower than the window would overflow its binning grid the same way.
! Runs AFTER the beam is built (the window length is the subject). setup_fel_elements
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
! Routine to set up chicane breaks (manual sec:phasing): a break whose elements bend the
! reference (sbends, patches) detours the BEAM while the RADIATION goes straight. The
! light's path is the CHORD between the flanking undulator faces, from ele%floor,
! never the reference arc that vec(5) is measured against. The arc-minus-chord
! delay is charged as whole-wavelength window rotations (Genesis's chicane
! semantics: "always autophasing") on the break's last element, which also takes
! the light-path drift correction. Absolute mode adds the delay's carrier phase in
! the walk. Only a closed bump keeps the light on the next undulator's axis.
! Anything else is refused by name, as is any geometry element under the
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
! Routine to set up the diagnostics (manual sec:stats): resolve the dump-at lists through
! Bmad's own lat_ele_locator (class::name syntax for free) and precompute the EXACT
! record and element-end counts by replaying the walk's skip rule, so the stats arrays
! are sized once, never grown. An entry matching nothing is refused by name.
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

! The record count REPLAYS the walk's skip rule, its window, its z arithmetic and
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
allocate (run%bdiag_arr(nslice), run%fpow_arr(nslice), run%fonax_arr(nslice))

end subroutine setup_diagnostics

end subroutine fel_setup_schedule


end module fel_setup_mod
