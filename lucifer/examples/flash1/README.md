# FLASH1 at 13.7 nm: a real machine at its published parameters

One command, Bmad only:

    ../../../production/bin/lucifer lucifer.in
    python ../plot_fel.py flash1.stats.h5

FLASH at DESY lased at 13.7 nm in 2006, the first free-electron laser to reach saturation at that wavelength, and reported its undulator, its beam and its performance. This example is that line: six planar segments, a Gaussian bunch on the published beam, and SASE from the quiet start through 27 m of undulator. Every other example in this directory runs the Aramis benchmark at 0.1 nm, so this one is also the third machine the numerical defaults are measured on ([SASE convergence](../../doc/startup-noise.md)).

The two sources are Ackermann et al., *Operation of a free-electron laser from the extreme ultraviolet to the water window*, Nature Photonics **1**, 336 (2007), [doi:10.1038/nphoton.2007.76](https://doi.org/10.1038/nphoton.2007.76), and Schreiber and Faatz, *The free-electron laser FLASH*, High Power Laser Science and Engineering **3**, e20 (2015), [doi:10.1017/hpl.2015.16](https://doi.org/10.1017/hpl.2015.16).

## What is published and what this file chose

| quantity | published | source | this file | note |
|---|---|---|---|---|
| undulator period | 2.73 cm | Ackermann, results | 0.0273 m | |
| peak field | 0.47 T | Ackermann, results | `b_max = 0.47` | |
| undulator parameter | K = 1.23, K_rms = 0.9 | Schreiber 3.4 | K = 1.19807, aw = 0.847162 | derived from the field, 2.7 percent under the published K |
| gap | 12 mm fixed | both | not modelled | a Bmad wiggler is described by its field |
| segments | six, 4.5 m each | Ackermann, results | six, 4.5045 m | 165 whole periods, since 4.5 m is 164.84 of them |
| intersection | doublet, monitors, wire scanners | Schreiber 3.4 | 0.8 m, `QF` and `QD` at k1 = 1.3 and -0.8 | length and strengths chosen here |
| beta function | about 10 m | Schreiber 3.4 | 10.0 m mean, 9.3 to 10.6 m | what the chosen doublet gives |
| beam energy | 700 MeV, spike near 680 MeV | Ackermann, results and Fig. 2 | 668.494 MeV | the resonance at 13.7 nm for K = 1.19807 |
| peak current | 2 to 2.5 kA | Ackermann, saturated output | 2.5 kA | the upper end |
| spike length | about 30 fs FWHM | Ackermann, saturated output | `sig_z = 3.8193e-6` m | a Gaussian of that width, 79.8 pC |
| normalized emittance | 1 to 1.5 mm mrad | Ackermann, saturated output | 1.5 µm both planes | the upper end |
| energy spread | below 0.1 percent | Schreiber 2 | `sig_pz = 6e-4` | argued below |
| Pierce parameter | 2.5 to 3 x 10^-3 | Ackermann, saturated output | 1.97 x 10^-3 | this beam, with the planar coupling |

Three numbers are this file's and are marked in the lattice header. The segment is a whole number of periods because the FEL step wants one. The intersection has a length because a drift must, and the papers give none. The doublet has strengths because the beta function of 10 m has to come from somewhere. A planar undulator focuses vertically only, and its own focusing gives 6.71 m vertically and nothing horizontally, so the doublet makes the horizontal beta function and raises the vertical one.

The published K and the published field disagree by 2.7 percent. K = 1.23 needs 0.4825 T at this period. The field is the one taken, and the disagreement then appears again in the energy: 13.7 nm needs 668.5 MeV at K = 1.19807 and 686 MeV at K = 1.23, against the machine's 700 MeV with its lasing spike near 680.

The slice energy spread. The papers give a bound on the beam and no slice value. Figure 2 of Ackermann puts a chirp of about 0.1 MeV per fs across the lasing spike. The FEL averages over one coherence time, measured there as 4.2 fs, so the chirp contributes 0.42 MeV to the spread the interaction sees, which is 6e-4 of the energy. That is the value here, and it is under the published bound of 0.1 percent. It matters: at 1e-4 this deck's field gain length is 1.6 m and at 1e-3 it is 2.6 m.

## What the run measures

Four seeds at the settings below. The spread is the SASE fluctuation of a 351-slice window.

| quantity | this deck | published | source |
|---|---|---|---|
| field gain length in the exponential regime | 2.05 +- 0.16 m | 2.5 +- 0.3 m | Ackermann Fig. 4a |
| end of the exponential regime | 20.9 +- 1.1 m | inside the 27 m line | Ackermann Fig. 4a |
| pulse energy there | 26.6 +- 5.7 µJ | 40 µJ characterized, 70 µJ demonstrated | Ackermann, results |
| pulse energy at the undulator exit | 50.3 +- 3.1 µJ | | |
| peak power in one slice at the exit | 2.8 to 5.0 GW | 10 GW | Ackermann, abstract |
| peak bunching factor | 0.244 +- 0.033 | | |

The gain lengths agree where their bands meet and the pulse energies agree. The simulation reaches the exponential regime sooner than the machine did, which is what an ideal beam does. This deck has no energy chirp along the bunch, no wakefield from the 12 mm gap, and no error in the trajectory or the matching, and each of those lengthens a real gain curve. The chirp is the largest of them. Ackermann reports a 10 fs radiation pulse out of a 30 fs spike, because the chirp takes the rest of the spike off resonance. This deck lases over the whole spike, so its pulse is longer than the machine's.

The gain length is read from the 4 m window of fastest growth in the pulse energy. Read instead over the whole rise from startup to saturation it is 2.9 m, since the curve is not one exponential. Ackermann fit the exponential part of Fig. 4a, which is the same choice.

## The numerical settings, and how each was chosen

Every value is measured on this machine. The Aramis column is what the other examples use at 0.1 nm on a beam five times smaller.

| setting | Aramis | here | measured |
|---|---|---|---|
| `grid_n_pts` | 255 | 129 | cells of sigma_x/6.8, the criterion of the convergence page. 65 points give 44.9 µJ against 46.3, and 257 raise the power outside the mode from 7 to 27 percent |
| `grid_half_width` | 2e-4 m | 1.0e-3 m | 9.3 rms beam sizes. At the same cell size 0.7 mm and 1.5 mm give 46.29 and 46.35 µJ, so the half width does not enter |
| `ds_step` | 0.045 m | 0.0819 m | three periods, 55 steps per segment. One period gives 48.8 µJ against 46.3, inside the seed spread, at 2.7 times the wall clock |
| `slicing%n_wavelength` | 3 | 12 | 3.4 slices per cooperation length of 40 wavelengths. Six wavelengths per slice give 48.5 µJ against 46.3, inside the seed spread, at twice the wall clock |
| `slicing%n_slice` | 96 | 351 | the bunch to four sigma either side, plus 990 wavelengths of slippage at each end |
| `beam_init%n_particle` | 2048 | 4096 | 1024 leaves 25 percent of the exit power outside the mode, 4096 leaves 7 percent, 16384 leaves 2 percent. The convergence page asks for under 10 |
| `beamlet_size` | 8 | 8 | at 4096 particles, 4 leaves 3 percent outside the mode and 16 leaves 11. Eight resolves three harmonics where four resolves one |
| `global%source_filter` | F | F | on at its derived 81 µrad edge it removes 162 times the wide-angle power and moves the pulse energy 0.8 percent, since 4096 particles has already converged. At 1024 particles it removes 219 times and reaches the same answer at two thirds of the wall clock |
| `resample%slice_width` | 0.01 | not read | the deck describes its bunch and the loader evaluates it per slice, so nothing is resampled |

Two of these do not carry over from Aramis. The convergence page estimates the wide-angle share in advance as 100/N_b times (1.57 µm/dx) squared, which gives 0.2 percent here where the measurement gives 7.9. The prefactor is the first machine's, as the page says, and it is 40 times low on this one, over cell sizes of 7.8 to 31 µm. The other is the cell size. The criterion sigma_x/7 holds, and it lands on 129 points at a 1 mm half width rather than 255 at 0.2 mm, so a deck copied from `../sase` would have run a grid 16 times too fine and reported a power that is mostly artifact.

## The beam and the window

The lattice carries the optics, so the bunch is generated matched to the periodic solution in the `beginning` statement, and there is no `gamma0` and no match transform. The current is derived from `bunch_charge` and `sig_z` and is never given. The window covers the bunch to four sigma either side with 990 wavelengths of slippage headroom at each end, because radiation made at the head of the bunch slips out of a shorter window before the exit. About a quarter of the slices therefore carry no charge, which is the price of holding the whole pulse.

Measured on this input at the `ran_seed` default of 12345: the exit power is 84.5 GW summed over the window and the pulse energy is 46.3 µJ. The run's own convergence report puts the power outside the split angle of 81 µrad at 0.078 times the power inside it, so 7 percent of that total is the wide-angle emission of the point beamlets and the rest is the mode.

Leaving `grid_n_pts` and `grid_half_width` out of the deck derives 127 points over a 962 µm half width, cells of 15.27 µm against the 15.63 µm here. That run gives 46.33 µJ against 46.35 and a ratio of 0.080 against 0.078, so the stated grid is the derived one.

Runs in ~20 s.
