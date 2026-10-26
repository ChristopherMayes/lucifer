# Longitudinal space charge

The example runs in four commands, with Bmad only:

    ../../build-release/bin/lucifer lucifer.toml
    ../../build-release/bin/lucifer lucifer_off.toml
    ../../build-release/bin/lucifer lucifer_short_range.toml
    ../../build-release/bin/lucifer lucifer_gate_off.toml
    python -m lucifer.plot space_charge.stats.h5

Space charge enters the pendulum equation as a per-particle `ez`, from per-slice radial-harmonic solves plus a whole-window long-range term (the manual's [space-charge section](../../docs/fel-physics.md)). The two terms answer different questions, so this example measures them separately.

The elements say where it acts. `sc_line.bmad` calls the shared lattice and sets Bmad's own attribute on the undulators, `wiggler::*[SPACE_CHARGE_METHOD] = slice`, and the input files turn on Bmad's master switch, `bmad.csr_and_space_charge_on`. The FEL slices are the bins the solve uses, so Bmad's `space_charge_com` `n_bin`, set in the lattice file plays no part in it. `space_charge.` in the input file holds the solver's numbers and nothing about whether it runs.

`lucifer.toml` and `lucifer_off.toml` isolate the long-range term on a cold dark beam: a 4 nm Gaussian bunch at 3 kA peak, no seed and no shot noise, so space charge is the only channel that can change any particle's energy. The control run is the proof of that. With space charge off, every slice's mean energy is unchanged to the last bit over 57 m.

The measured quantities of the four input files are in the table below, the head-to-tail span of the slice mean energies at the exit among them, with the control run's exact zero beside it. The long-range term gives all of the span, and the two short-range harmonics the input file sets give none of it.

The beam-energy panel of the plot is the whole result: a symmetric fan opening from zero to half the span either side, with the mean and the energy spread both flat on zero through all 57 m. The bunch pushes itself apart, the tail losing energy and the head gaining it, and the mean is preserved because the force is internal. The power panels show femtowatt numbers, which are the dark run's numerical floor and not light. A 4 nm bunch is 13 attoseconds long, which is where this matters. Lengthening it ten times at the same 3 kA peak, by scaling `sig_z`, `bunch_charge` and `window_length` together, drops the span to about a third. The effect weakens with bunch length and not in proportion to it. A 4 nm bunch is far shorter than the beam's own 20 micrometer transverse size, and the simple line-charge scaling does not hold there.

The chirp matters at this size. Relative to the beam energy the span is more than twenty times the bandwidth that the gain length of `../steady_state` implies, so a bunch this short does not lase without something to compensate it.

The short-range harmonics measure exactly zero here, which is correct. They act on the density structure inside a slice, and a quiet-start dark beam has none, since the initial bunching sits at the quiet start's floor. `lucifer_short_range.toml` gives them something to act on, by matching the `../steady_state` input file parameter for parameter and adding `nz = 4, nphi = 1`. Against that run as its control, whose exit power and bunching its own page carries, the short-range input file's exit power is higher by about half a percent. That rise comes past saturation, from the short-range term acting on a beam the FEL has bunched.

Two configurations are refused. `slice` with neither term configured, `nz = 0` and `longrange = false`, is refused, since the solve would run at its full cost and return an exact zero. `fft_3d` on an FEL element is refused, since that solver wants a three-dimensional grid this walk does not build.

The master switch off while the elements still ask for `slice` is not refused, because that is how a lattice is run without space charge and without editing it. `lucifer_gate_off.toml` is that case: the same lattice, the same solver numbers, the switch false. The run prints a note saying so and tracks without space charge, measured as exactly the control run's zero.

The run reports the split between the mode and the wide-angle emission of the point beamlets. That emission comes from macroparticles that occupy one grid point each, and no real beam radiates it ([SASE convergence](../../docs/startup-noise.md)). At the exit the power outside the split angle is a fraction of a percent of the power inside it on the short-range input file, so its exit power is the mode, as is the control's. With the filter off both reported nearly half of their power as that emission. The long-range input files start dark, and what the filter leaves of their femtowatt fields is a few percent outside the mode, which is why the comparison there is made on the beam.

Each run takes about seven seconds.
