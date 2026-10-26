"""The input file as typed models, generated from the schema. Do not edit.

Written by python -m lucifer.generate from schema.json, schema version 1. A dataclass per table of the input file, and Input with
the lattice and one of each. The defaults are the program's, read from the schema,
and nothing here checks a bound: the program refuses a value outside one. The
field named field is the [field] table, so dataclasses.field is imported as _field.
"""
from __future__ import annotations

from dataclasses import dataclass, field as _field

from lucifer.inputops import InputOps

__all__ = ["Input", "Run", "Window", "Beam", "Resample", "Field", "Wake", "SpaceCharge", "Bmad"]


@dataclass
class Run:
    """The [run] table of the input file."""
    out_root: str = 'fel_track'
    """Root of every output file name."""
    interlude_model: str = 'bmad'
    """Interlude transport: bmad (the seam) or genesis (transcribed)."""
    transport_model: str = 'bmad'
    """Transverse transport inside averaged FEL elements: bmad, Bmad's own kernel, or genesis, the transcribed maps for validation only."""
    unaveraged_steps_per_period: int = 20
    """Integration steps per undulator period in unaveraged elements. Below 10 the push does not converge. Bounds: from 10."""
    unaveraged_ramp_periods: int = 2
    """Periods the field ramps over at each end of an unaveraged element. 0 keeps the default of 2, -1 is a hard edge. Bounds: from -1."""
    source_model: str = 'deposit'
    """The FEL source model: deposit, the per-particle scatter, or coherent, the coherent-Gaussian source."""
    source_filter: bool = True
    """The angular filter on the source term, transcribed from Genesis4. Off to compare against a code without one."""
    source_filter_angle: float = 0.0
    """The filter's edge as an angle. 0 derives it from the beam and the gain. [rad] Bounds: from 0.0."""
    source_filter_tolerance: float = 0.75
    """The largest fraction of source intensity the filter may remove inside the protected angle. Bounds: above 0.0 below 1.0."""
    source_filter_xcut: float = 0.0
    """The edge in Genesis4's own unit, the shifted grid index over ngrid. Validation only. Bounds: from 0.0."""
    source_filter_ycut: float = 0.0
    """The edge in Genesis4's own unit, the shifted grid index over ngrid. Validation only. Bounds: from 0.0."""
    source_filter_width: float = 0.05
    """The sigmoid's width as a fraction of its edge. Small is a sharp edge. Bounds: above 0.0."""
    source_filter_mutate: bool = False
    """Self-test of the filter check: apply the sigmoid to the propagated field instead of the source."""
    track_start: str = ''
    """Element locator opening the walk. Blank is the first element."""
    track_end: str = ''
    """Element locator closing the walk. Blank is the last element."""
    continuation: bool = False
    """beam.file and field.file are a checkpoint pair the run continues from, not an initialization."""
    comb_ds_save: float = 0.0
    """Minimum z advance between per-record stats rows. 0 is every record, negative is none. [m]"""
    interlude_ds_step: float = 0.0
    """The length an interlude element is cut into, for stats rows and frames along it. 0 tracks it whole. [m]"""
    dump_beam_at: list[str] = _field(default_factory=lambda: [])
    """Element locators for mid-run beam dumps."""
    dump_field_at: list[str] = _field(default_factory=lambda: [])
    """Element locators for mid-run field dumps."""
    dump_at_comb: bool = False
    """Write the beam and the field at every comb position as a frame series."""
    dump_slice_first: int = 1
    """The first slice a frame carries, in window order. Bounds: from 1."""
    dump_slice_last: int = -1
    """The last slice a frame carries. -1 is the last slice of the window. Bounds: from -1."""
    dump_reduced: bool = False
    """Write the field's own reductions per record into the stats file."""
    ran_seed: int = 12345
    """The one random seed: generation, import and noise."""
    write_diag: bool = False
    """Write the Genesis-comparison text diag file, which is large."""
    write_initial: bool = False
    """Dump the initial state before tracking."""
    load_only: bool = False
    """Build the initial state, dump it, stop."""
    keep_escaped_field: bool = True
    """Keep the light slippage carries out of the window, as the escaped bank and the whole pulse at the exit. On by default, since the light is what an FEL makes. A run whose window holds all its light writes neither file."""
    migrate: bool = False
    """Move particles between slices when their phase leaves the slice window."""
    migrate_check: bool = False
    """The migration's bunching-invariance instrument."""
    reference_run: bool = False
    """No FEL interaction: Bmad tracks everything."""
    record_environment: bool = False
    """Record the login name and the working directory in meta/. Off, since a stats file travels."""
    fp32_check: str = 'off'
    """The FP32 lockstep check: off, lockstep or freerun. For checks only."""
    fp32_mutate: bool = False
    """The FP32 check's self-test: coarsen the residual."""
    device: str = 'off'
    """The device backend: off, or metal on a build that carries it."""
    device_timing: bool = False
    """Per-pass device timing, one compute encoder a pass. Slows the run."""
    device_dep_mutate: int = 0
    """Shift the device deposit's fixed-point scale, in bits. A check instrument."""


@dataclass
class Window:
    """The [window] table of the input file."""
    window_length: float = 0.0
    """The time window. 0 derives it from the bunch and the slippage the line applies, or is one slice. [m] Bounds: from 0.0."""
    n_wavelength: int = 0
    """Wavelengths per slice, Genesis4's sample. 0 derives it from the gain. Bounds: from 0."""
    n_slice: int = 0
    """The window as a slice count. 0 takes it from window_length. Both set is refused. Bounds: from 0."""
    current: float = 0.0
    """A flat current. 0 derives the current from the bunch's z structure. [A] Bounds: from 0.0."""


@dataclass
class Resample:
    """The [beam.resample] table of the input file."""
    slice_width: float = 0.01
    """Sampling window over bunch length, Genesis4's slicewidth. Bounds: above 0.0."""
    n_particle_per_slice: int = 8192
    """Macroparticles per slice on the validation route. Bounds: from 1."""
    beamlet_size: int = 4
    """Beamlet size on the validation route. Bounds: from 1."""
    n_slice: int = 0
    """Slice count. 0 derives it from the bunch length and the spacing. Bounds: from 0."""
    use_beam_init: bool = False
    """Validation route: resample a bunch Bmad generated from [beam] instead of the analytic loader."""


@dataclass
class Beam:
    """The [beam] table of the input file."""
    file: str = ''
    """An openPMD particle dump to start from, already in slices. Blank loads a bunch instead."""
    position_file: str = ''
    """An openPMD-beamphysics bunch to load, read whole and then sliced by load_mode."""
    n_particle: int = 0
    """Macroparticles per slice after the phase copies, a positive multiple of beamlet_size. In the keep mode, the particles of the generated bunch. Bounds: from 1."""
    bunch_charge: float = 0.0
    """Total bunch charge. The current is derived from it and never input. [C] Bounds: from 0.0."""
    sig_z: float = 0.0
    """Bunch length, read with z_distribution. 0 is the steady state of one slice. [m] Bounds: from 0.0."""
    sig_pz: float = 0.0
    """Fractional momentum spread dP/P0. Bounds: from 0.0."""
    a_norm_emit: float = 0.0
    """Normalized a-mode emittance. [m rad] Bounds: from 0.0."""
    b_norm_emit: float = 0.0
    """Normalized b-mode emittance. [m rad] Bounds: from 0.0."""
    z_distribution: str = 'RAN_GAUSS'
    """The current profile: ran_gauss, a Gaussian of sig_z, or grid, flat from z_min to z_max."""
    z_min: float = 0.0
    """The flat profile's lower z bound, read under grid only. [m]"""
    z_max: float = 0.0
    """The flat profile's upper z bound, read under grid only. [m]"""
    load_mode: str = 'sample'
    """How the bunch becomes slices: sample, the same count in every slice with weights from the slice charge, or keep, every particle kept with the counts following the charge."""
    quiet_start: bool = True
    """Load beamlets of beamlet_size phase copies, quiet below that harmonic."""
    beamlet_size: int = 8
    """Copies per beamlet, and the harmonics the noise resolves. Bounds: from 1."""
    shot_noise: bool = False
    """Impose the physical noise on the load's groups. Time-dependent windows only."""
    write_genesis_dist: str = ''
    """Write the bunch as a Genesis4 importdistribution file."""
    write_openpmd_file: str = ''
    """Write the bunch as openPMD-beamphysics."""
    split_weights: bool = False
    """Check instrument: coincident copies carrying 1/3 and 2/3 of each weight after loading."""
    swap_beam_xy: bool = False
    """Check instrument: swap (x, px) with (y, py) after generation."""
    gen_test_weights: bool = False
    """Check instrument: alternate beamlet weights of 0.25 and 1.75, charge preserving."""
    resample_split_weights: bool = False
    """Check instrument: the split-weight copies applied to the bunch before it is sliced."""
    resample: Resample = _field(default_factory=Resample)
    """The [beam.resample] table, the resampler."""


@dataclass
class Field:
    """The [field] table of the input file."""
    file: list[str] = _field(default_factory=lambda: [])
    """The openPMD wavefront file to start from, or one per harmonic member. Blank is a seed."""
    lambda0: float = 0.0
    """The resonant wavelength. Required for generation. [m] Bounds: from 0.0."""
    grid_n_pts: int = 0
    """Transverse grid points per side. 0 derives it from the beam. Bounds: from 0."""
    grid_half_width: float = 0.0
    """Transverse grid half width. 0 derives it from the beam. [m] Bounds: from 0.0."""
    seed_power: float = 0.0
    """Gaussian seed power. 0 is a dark start. [W] Bounds: from 0.0."""
    seed_waist_size: float = 0.0
    """The seed's intensity 1/e^2 radius. [m] Bounds: from 0.0."""
    seed_polarization: str = 'x'
    """The seed's polarization, x or y."""
    harmonics: list[int] = _field(default_factory=lambda: [1])
    """The field set: 1 first, then further harmonic numbers in increasing order."""


@dataclass
class Wake:
    """The [wake] table of the input file."""
    on: bool = False
    """The chamber wake."""
    model: str = 'genesis'
    """Which chamber-wake implementation runs: genesis, the transcribed solver."""
    loss: float = 0.0
    """External loss. [eV/m]"""
    radius: float = 0.0025
    """Chamber radius, or the half gap of a flat chamber. [m] Bounds: above 0.0."""
    conductivity: float = 0.0
    """DC conductivity. 0 is no resistive wake. [1/(Ohm m)] Bounds: from 0.0."""
    relaxation: float = 0.0
    """AC relaxation distance c tau. [m] Bounds: from 0.0."""
    roundpipe: bool = True
    """A round chamber. False is flat, parallel plates."""
    material: str = ''
    """CU or AL, a shortcut for the conductivity and relaxation."""
    gap: float = 0.0
    """Undulator gap. 0 is no geometric wake. [m] Bounds: from 0.0."""
    lgap: float = 1.0
    """Period of the gaps. [m] Bounds: above 0.0."""
    hrough: float = 0.0
    """Roughness amplitude. 0 is no roughness wake. [m] Bounds: from 0.0."""
    lrough: float = 1.0
    """Roughness period. [m] Bounds: above 0.0."""
    write_kernels: str = ''
    """Check instrument: write the transcribed kernels to this file."""


@dataclass
class SpaceCharge:
    """The [space_charge] table of the input file."""
    model: str = 'genesis'
    """Which space-charge implementation runs: genesis, the transcribed solver."""
    rmax: float = 0.0
    """Radial grid extent scale. It grows adaptively as Genesis4's does. [m] Bounds: from 0.0."""
    ngrid: int = 100
    """Radial grid points. Bounds: from 1."""
    nz: int = 0
    """Longitudinal harmonics. 0 disables the short-range solve. Bounds: from 0."""
    nphi: int = 0
    """Azimuthal modes, -nphi to nphi. Bounds: from 0."""
    longrange: bool = False
    """The whole-window long-range term."""


@dataclass
class Bmad:
    """The [bmad] table of the input file."""
    radiation_damping_on: bool = False
    """Bmad's radiation damping in the interludes."""
    radiation_fluctuations_on: bool = False
    """Bmad's radiation fluctuations in the interludes."""
    csr_and_space_charge_on: bool = False
    """Bmad's master switch for CSR and space charge."""
    sr_wakes_on: bool = True
    """Bmad's short-range element wakes."""
    lr_wakes_on: bool = True
    """Bmad's long-range element wakes."""
    default_ds_step: float = 0.2
    """Bmad's default tracking step. [m] Bounds: from 0.0."""
    space_charge_mesh_size: list[int] = _field(default_factory=lambda: [32, 32, 64])
    """Bmad's space-charge mesh, three counts."""
    sc_min_in_bin: int = 10
    """Bmad's minimum particles per space-charge bin. Bounds: from 1."""
    rel_tol_tracking: float = 1e-08
    """Bmad's relative tracking tolerance. Bounds: from 0.0."""
    particle_sigma_cutoff: float = -1.0
    """Bmad's particle sigma cutoff for the space-charge mesh."""


@dataclass(repr=False)
class Input(InputOps):
    """One input file: the lattice and the tables."""
    lattice: str = ''
    """The Bmad lattice file. The FEL segments are its wigglers with fel_method set."""
    run: Run = _field(default_factory=Run)
    """The [run] table."""
    window: Window = _field(default_factory=Window)
    """The [window] table."""
    beam: Beam = _field(default_factory=Beam)
    """The [beam] table."""
    field: Field = _field(default_factory=Field)
    """The [field] table."""
    wake: Wake = _field(default_factory=Wake)
    """The [wake] table."""
    space_charge: SpaceCharge = _field(default_factory=SpaceCharge)
    """The [space_charge] table."""
    bmad: Bmad = _field(default_factory=Bmad)
    """The [bmad] table."""


TABLES = (("run", Run), ("window", Window), ("beam", Beam), ("beam.resample", Resample), ("field", Field), ("wake", Wake), ("space_charge", SpaceCharge), ("bmad", Bmad))
"""The tables in the input file's order, with their classes."""
