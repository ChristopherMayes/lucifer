!+
! Module fel_setup_mod
!
! The lattice walk-in of the FEL tracker: hooks and lattice attributes, the parse,
! element recognition with its refusals, the harmonics/field-set validation, the
! collective configuration, the slippage/autophasing schedule, geometry breaks and
! the diagnostics setup. Split as fel_setup_lattice (needs only the inputs) and
! fel_setup_schedule (needs the built starting state). Library contract: errors
! return through err_flag; nothing here stops. The print lines are unchanged from
! when this code lived in the driver.
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
!+
! Subroutine fel_setup_lattice (run, err_flag)
!
! Everything that turns the parsed inputs into a recognized lattice: the tracking
! hooks and FEL lattice attributes (registered once per process), the parse itself,
! element recognition and its refusals (setup_fel_elements), the harmonics-request
! validation and the field-set allocation. Errors return through err_flag; nothing
! here stops (the library contract).
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
character(8) wavefront_format
character(1) seed_polarization
logical migrate, reference_run, err
integer :: harmonics(9)
integer n_harm, ih

!

err_flag = .false.
lat => run%lat
lat_file = run%lat_file
interlude_model = run%global%interlude_model
wavefront_format = run%global%wavefront_format
seed_polarization = run%winit%seed_polarization
migrate = run%global%migrate
reference_run = run%global%reference_run
harmonics = run%winit%harmonics

! (The FEL tracking mode is a lattice attribute; its refusals live in
! setup_fel_elements, and the unaveraged-vs-collective refusal follows setup, once
! any_unavg is known.)

if (interlude_model /= 'bmad' .and. interlude_model /= 'genesis') then
  print '(a)', 'fel_track_test: interlude_model must be "bmad" or "genesis", got: ' // trim(interlude_model)
  err_flag = .true.;  return
endif

! (bmad_com%radiation_damping_on / %radiation_fluctuations_on come straight from the
! &fel_params namelist -- Bmad's own switches, exposed directly as Tao exposes them.)

if (bmad_com%radiation_fluctuations_on .and. migrate) then
  print '(a)', 'fel_track_test: radiation_fluctuations draws one kick per BEAMLET (the quiet start'
  print '(a)', '  cancels per beamlet), and slice migration scrambles beamlet grouping. Pick one.'
  err_flag = .true.;  return
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

! The FEL mode and unaveraged parameters live on the LATTICE (Stage A of the CUSTOM
! retirement; manual sec:element), registered program-side so no lattice declares
! them. The same slot index serves wigglers and undulators.

if (.not. fel_attributes_registered) then
  call set_custom_attribute_name ('WIGGLER::FEL_TRACKING', err, 1)
  if (.not. err) call set_custom_attribute_name ('UNDULATOR::FEL_TRACKING', err, 1)
  if (.not. err) call set_custom_attribute_name ('WIGGLER::FEL_STEPS_PER_PERIOD', err, 2)
  if (.not. err) call set_custom_attribute_name ('UNDULATOR::FEL_STEPS_PER_PERIOD', err, 2)
  if (.not. err) call set_custom_attribute_name ('WIGGLER::FEL_RAMP_PERIODS', err, 3)
  if (.not. err) call set_custom_attribute_name ('UNDULATOR::FEL_RAMP_PERIODS', err, 3)
  if (err) then
    print '(a)', 'fel_track_test: could not register the FEL lattice attributes.'
    err_flag = .true.;  return
  endif
  fel_attributes_registered = .true.
endif

! err_flag matters: bmad_parser reports attribute errors (e.g. a wake on an element
! type that cannot carry one) and RETURNS -- without the check the run continues on a
! partial lattice, found when an example's drift wakes silently never attached.

call bmad_parser (lat_file, lat, err_flag = err)
if (err) then
  print '(a)', 'fel_track_test: lattice parse errors (above); refusing to run on a partial lattice.'
  err_flag = .true.;  return
endif
branch => lat%branch(0)

run%gamma0 = branch%ele(0)%value(e_tot$) / m_electron
print '(a, f0.6, a)', 'fel_track_test: gamma0 = ', run%gamma0, ' (from the lattice e_tot).'

call setup_fel_elements ()   ! Needs only the parsed lattice; sets run%two_pol for the
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
  print '(a)', 'fel_track_test: harmonics(1) must be 1 -- the FUNDAMENTAL anchors the phase,'
  print '(a)', '  the phi0 advance and the slippage schedule; harmonic fields ride on it.'
  err_flag = .true.;  return
endif
do ih = 2, 9
  if (harmonics(ih) == 0) then
    if (any(harmonics(ih:) /= 0)) then
      print '(a)', 'fel_track_test: harmonics must be a gap-free increasing list (0 padding at the end).'
      err_flag = .true.;  return
    endif
    exit
  endif
  if (harmonics(ih) <= harmonics(ih-1)) then
    print '(a)', 'fel_track_test: harmonics must be strictly increasing (no duplicates).'
    err_flag = .true.;  return
  endif
enddo
if (n_harm > 1 .and. any(fel_mode == fel_unaveraged$ .and. is_fel)) then
  print '(a)', 'fel_track_test: harmonic FIELDS with an UNAVERAGED element are not implemented'
  print '(a)', '  (the unaveraged mode carries the fundamental envelope only; its harmonic'
  print '(a)', '  couplings are validated through the particle spectra, not a carried field).'
  err_flag = .true.;  return
endif
if (n_harm > 1 .and. run%two_pol) then
  print '(a)', 'fel_track_test: harmonic fields with TWO LIVE POLARIZATIONS are not validated'
  print '(a)', '  together yet; run one or the other.'
  err_flag = .true.;  return
endif
select case (wavefront_format)
case ('genesis', 'openpmd', 'both')
case default
  print '(a)', 'fel_track_test: wavefront_format must be genesis, openpmd or both, not "' // &
                trim(wavefront_format) // '".'
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

subroutine setup_fel_elements ()

! Recognize FEL segments -- wiggler/undulator elements with tracking_method = custom,
! Bmad's own semantics for program-supplied tracking, which this driver is -- and derive
! their FEL parameters from lattice attributes, enforcing the brief's 7.5 assertions.
! The stored k1x/k1y wiggler attributes are deliberately NOT read (7.5: their helical
! sign disagrees with the tracking locals; nothing here cross-uses them).

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
  ! (Bmad's own kernel is the default); 1 = unaveraged; -1 = averaged with the
  ! transcribed-Genesis maps (validation-internal: the Genesis tiers require
  ! transcription-level transport; no production lattice writes it).

  rv = value_of_attribute(w, 'FEL_TRACKING', err_a)
  if (err_a) then
    err_flag = .true.;  return
  endif
  fel_mode(je) = nint(rv)
  if (abs(rv - fel_mode(je)) > 1e-9_rp .or. fel_mode(je) < fel_transcribed$ .or. &
      fel_mode(je) > fel_unaveraged$) then
    print '(a)', 'fel_track_test: fel_tracking must be -1 (transcribed maps, validation-internal),'
    print '(a)', '  0/unset (averaged, bmad_standard kernel maps) or 1 (unaveraged), at element: ' // trim(w%name)
    err_flag = .true.;  return
  endif

  rv = value_of_attribute(w, 'FEL_STEPS_PER_PERIOD', err_a)
  if (err_a) then
    err_flag = .true.;  return
  endif
  fel_spp(je) = nint(rv)
  if (fel_spp(je) == 0) fel_spp(je) = 20
  if (fel_spp(je) < 10) then
    print '(a)', 'fel_track_test: fel_steps_per_period is below the floor of 10 (MINERVA''s envelope),'
    print '(a)', '  at element: ' // trim(w%name)
    err_flag = .true.;  return
  endif

  ! fel_ramp_periods: an attribute's unset value is 0, and a silent hard edge would
  ! reintroduce the K/gamma handoff hazard by omission -- so unset/0 means the default
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
    print '(a)', 'fel_track_test: fel_ramp_periods must be positive, 0/unset (default 2), or the'
    print '(a)', '  hard-edge test sentinel -1, at element: ' // trim(w%name)
    err_flag = .true.;  return
  else
    fel_ramp(je) = rv
  endif

  ! The 7.5 assertions live in fel_assert_wiggler_sane, ONE authority called from the
  ! track1/mat6 hooks (where they fire first, during the parse) and again here. Keeping
  ! a second copy inline was tried and rejected: redundant assertions mask the removal
  ! of either copy, which defeats mutation testing of the refusal checks.

  call fel_assert_wiggler_sane (w)

  is_fel(je) = .true.
  kw = twopi / w%value(l_period$)

  ! aw (rms, Genesis's convention) from the peak field:
  ! K = c*b_max/(k_u * m_e c^2), exactly and independent of the reference energy;
  ! helical aw = K, planar aw = K/sqrt(2). Focusing split: Genesis's defaults by
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

  ! Tilt: the wiggle-plane rotation (planar only -- a tilted helical is a no-op that
  ! reads as confusion, refused; the transcribed-Genesis maps know no tilt, refused).
  ! The polarization 2-vector on (Ex, Ey): planar (cos t, sin t); helical (1,-i)/sqrt2.

  und_of(je)%tilt = w%value(tilt_tot$)
  if (und_of(je)%tilt /= 0) then
    if (und_of(je)%helical) then
      print '(2a)', 'fel_track_test: tilt on a HELICAL FEL element is a rotation of a ', &
                    'circularly symmetric field -- a no-op that reads as a mistake: ' // trim(w%name)
      err_flag = .true.;  return
    endif
    if (fel_mode(je) == fel_transcribed$) then
      print '(2a)', 'fel_track_test: the transcribed-Genesis maps (fel_tracking = -1) know ', &
                    'no tilt (Genesis has none); use the default maps on: ' // trim(w%name)
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
  print '(a)', 'fel_track_test: seed_polarization must be "x" or "y".'
  err_flag = .true.;  return
endif

if (.not. any(is_fel) .and. .not. reference_run) then
  print '(a)', 'fel_track_test: the lattice has no FEL elements (wiggler/undulator with tracking_method = custom).'
  print '(a)', '  Set reference_run = T for a deliberate no-FEL run (Bmad tracks everything).'
  err_flag = .true.;  return
endif



end subroutine setup_fel_elements

end subroutine fel_setup_lattice

!------------------------------------------------------------------------------
!+
! Subroutine fel_setup_schedule (run, err_flag)
!
! Everything that needs BOTH the lattice and the built starting state: the collective
! configuration and wake kernels, the wake-window and unaveraged-collective refusals,
! the slippage/autophasing schedule, the geometry breaks (chicanes), and the
! diagnostics setup (dump locators; the exact record/element-end counts). Errors
! return through err_flag; nothing here stops.
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
character(60) :: dump_beam_at(40), dump_field_at(40)
logical wake_on, sc_longrange, timerun, two_pol, err
real(rp) sc_rmax, wake_loss, wake_radius, wake_conductivity, wake_relaxation
real(rp) wake_gap, wake_lgap, wake_hrough, wake_lrough
logical wake_roundpipe
character(8) wake_material
integer sc_ngrid, sc_nz, sc_nphi
real(rp) gamma0_ref, phase_rate, Lz
integer nslice, n_harm, ie, ih, prev_ie, nrec_stats, nend_stats
integer :: harmonics(9)

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
  if (err) then
    err_flag = .true.;  return
  endif
  call fel_wake_update (coll%wake, fbeam)
  open (newunit = run%iu_wake, file = trim(out_root) // '.wake.txt', action = 'write')
  call fel_write_wake_block (run, 0.0_rp)

  ! The single-particle kernels, for the deliverable-11 cross-validation: the same
  ! physical wake fed to Bmad's z_long machinery must reproduce these kernels'
  ! convolution. NOTE the s = 0 entries carry the Bane self-slice half factor
  ! (fel_wake_init halves them); a plain W(z) table wants the unhalved value.

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
      print '(a)', 'fel_track_test: wrote wake kernels: ' // trim(run%wake_init%write_kernels)
    end block
  endif
endif

! with slippage active; one slice is the steady state and fel_apply_slippage is a no-op.

timerun = (nslice > 1)
do ih = 1, n_harm       ! Same values for every field: one window, lockstep rotation
                        ! in fundamental-wavelength units (Genesis's one Control::sample).
  ffield(ih)%slip%timerun = timerun
  ffield(ih)%slip%sample = fbeam%slice_spacing / fbeam%wavelength
enddo
run%gamma0_ref = fel_gamma0(fbeam)
gamma0_ref = run%gamma0_ref

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
! (setup_fel_elements itself runs right after the parse, before the beam and field are
! built: two_pol -- does any element tilt? -- must be known when the field is made.)

call check_wake_window ()
if (err_flag) return
run%any_unavg = any(fel_mode == fel_unaveraged$ .and. is_fel)

! The unaveraged mode is a verification mode (fel-physics.tex sec:unaveraged): the
! collective terms are not wired into its step, and a mixed line would apply them in
! some segments and silently drop them in others. Refuse by name.

if (run%any_unavg .and. (wake_on .or. sc_nz >= 1 .or. sc_longrange)) then
  print '(a)', 'fel_track_test: wakes/space charge are NOT wired into the unaveraged mode'
  print '(a)', '  (a verification mode; see fel-physics.tex sec:unaveraged). Turn them off.'
  err_flag = .true.;  return
endif


! The rest of the schedule: drift autophasing. Interludes accumulate Lz; the last
! interlude before each undulator gets floor(Lz/(2*gamma0^2*lambda)) + 1 wavelengths
! (Lattice.cpp:171-174, guarded there by Lz > 0). The end-of-lattice fixup
! (Lattice.cpp:191-193) is UNGUARDED in Genesis: the last element always gets
! floor(Lz/(2*gamma0^2*lambda)) + 1, which is +1 even with no trailing interlude at all
! ("autophasing is applied in case for [a] second, succeeding run"). Transcribed as is --
! omitting that +1 leaves the field record one rotation short at the very end, found the
! hard way against the single-segment time-dependent run. (Citations kept AT THE LINES:
! this quirk's exactness matters here, at the call site; manual sec:slippage.)

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
      print '(2a)', 'fel_track_test: the z_offset of FEL element ', trim(ele%name)
      print '(a, es10.2, a, es10.2, a)', '   (', fel_zoff(ie), ' m) exceeds its upstream break (', Lz, ' m).'
      err_flag = .true.;  return
    endif
    if (fel_zoff(ie) /= 0 .and. prev_ie == 0) then
      print '(2a)', 'fel_track_test: the z_offset knob on the FIRST element ', trim(ele%name)
      print '(a)', '   has no upstream break to displace into; give the lattice a leading break.'
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

! The tracking window (global%track_start/track_end, Tao's names; Genesis zstop
! parity). Resolved through Bmad's own locator; the SCHEDULE above was built on the
! full lattice, so a windowed run composes exactly with the full one -- the walk
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
  print '(a, i0, a, i0, a)', 'fel_track_test: track_start (element ', run%i_start, &
        ') is past track_end (element ', run%i_end, ').'
  err_flag = .true.;  return
endif

call setup_diagnostics ()
if (err_flag) return

!------------------------------------------------------------------------------
contains

!------------------------------------------------------------------------------
! Resolve one tracking-window locator to a single tracked-element index, refused by
! name when it matches nothing or more than one element.

subroutine resolve_window_ele (locator, which, ix)

type (ele_pointer_struct), allocatable :: eles(:)
character(*) locator, which
integer ix, n_loc
logical lerr

call lat_ele_locator (locator, lat, eles, n_loc, lerr)
if (lerr .or. n_loc == 0) then
  print '(4a)', 'fel_track_test: ', which, ' matches no element: ', trim(locator)
  err_flag = .true.;  return
endif
if (n_loc > 1) then
  print '(4a)', 'fel_track_test: ', which, ' matches more than one element: ', trim(locator)
  err_flag = .true.;  return
endif
if (eles(1)%ele%ix_ele < 1 .or. eles(1)%ele%ix_ele > branch%n_ele_track) then
  print '(4a)', 'fel_track_test: ', which, ' is not a tracked element: ', trim(locator)
  err_flag = .true.;  return
endif
ix = eles(1)%ele%ix_ele

end subroutine resolve_window_ele

!------------------------------------------------------------------------------
! concatenate into one bunch in global window coordinates and Bmad's wake machinery
! applies unmodified. What is checked here, by name: lr (multi-bunch) wakes are not
! supported; a pseudomode wake whose z_max is shorter than the window would have Bmad
! kill the bunch mid-run; a z_long table narrower than the window would overflow its
! binning grid the same way. Runs AFTER the beam is built (the window length is the
! subject); setup_fel_elements runs before it (two_pol must precede the field).

subroutine check_wake_window ()

type (ele_struct), pointer :: w
integer je

do je = 1, branch%n_ele_track
  w => pointer_to_wake_ele(branch%ele(je))
  if (.not. associated(w)) cycle
  if (allocated(w%wake%lr%mode)) then
    if (size(w%wake%lr%mode) > 0) then
      print '(2a)', 'fel_track_test: lr (multi-bunch) wakes are not supported; ', &
                    'remove them from: ' // trim(w%name)
      err_flag = .true.;  return
    endif
  endif
  if (w%wake%sr%z_max > 0 .and. size(fbeam%slice) * fbeam%slice_spacing > w%wake%sr%z_max) then
    print '(2a)', 'fel_track_test: the time window is longer than this element''s sr wake ', &
                  'z_max can handle: ' // trim(w%name)
    err_flag = .true.;  return
  endif
  if (w%wake%sr%z_long%dz > 0 .and. &
      size(fbeam%slice) * fbeam%slice_spacing > w%wake%sr%z_long%z0) then
    print '(2a)', 'fel_track_test: the time window is longer than this element''s z_long ', &
                  'wake table extent z0: ' // trim(w%name)
    err_flag = .true.;  return
  endif
enddo

end subroutine check_wake_window

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Element sr wakes act across the WHOLE window (deliverable 11): all slices
! Chicane breaks (manual sec:phasing): a break whose elements bend the reference
! (sbends, patches) detours the BEAM while the RADIATION goes straight -- so the
! light's path is the CHORD between the flanking undulator faces, from ele%floor,
! never the reference arc that vec(5) is measured against. The arc-minus-chord
! delay is charged as whole-wavelength window rotations (Genesis's chicane
! semantics: "always autophasing") on the break's last element, which also takes
! the light-path drift correction; absolute mode adds the delay's carrier phase in
! the walk. Only a CLOSED BUMP keeps the light on the next undulator's axis --
! anything else is refused by name, as is any geometry element under the
! genesis-model interludes (Genesis's drift/quad set cannot represent it).

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
        print '(3a)', 'fel_track_test: geometry element ', trim(ele%name), ' (a bend or patch) inside a'
        print '(a)', '   GENESIS-MODEL interlude: Genesis''s drift/quad set cannot represent it.'
        print '(a)', '   Use interlude_model = "bmad" (the seam tracks it exactly).'
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
! One geometry break: the chord, the closed-bump refusal, the delay's rotations and
! the light-path correction (see setup_break_geometry's header).

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
  print '(3a)', 'fel_track_test: the break ending at ', trim(branch%ele(iblast)%name), ' bends the reference'
  print '(a)', '   and is NOT A CLOSED BUMP: the radiation would leave the next undulator''s axis.'
  print '(a)', '   Only closed-bump chicanes are modeled; close the geometry or straighten the line.'
  err_flag = .true.;  return
endif

delay_geo = arc - chord
ele_slip(iblast) = ele_slip(iblast) + floor(delay_geo / wf%wavelength)
light_corr(iblast) = delay_geo

end subroutine close_geometry_break

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Diagnostics (manual sec:stats). setup_diagnostics resolves the dump-at lists through
! Bmad's own lat_ele_locator (class::name syntax for free; an entry matching nothing is
! refused by name) and precomputes the EXACT record and element-end counts by replaying
! the walk's skip rule -- the stats arrays are sized once, never grown.

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
    print '(2a)', 'fel_track_test: dump_beam_at entry matches no element: ', trim(dump_beam_at(i))
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
    print '(2a)', 'fel_track_test: dump_field_at entry matches no element: ', trim(dump_field_at(i))
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
