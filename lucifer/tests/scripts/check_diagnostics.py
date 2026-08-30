#!/usr/bin/env python3
"""
Diagnostic-output checks (manual sec:stats): the stats file, dump-at elements, and the
escaped-field bank, held by cross-identities rather than reference files --

  1. bunch_params reconstruction: bunch_params_from_stats.py applied to the per-record
     sufficient statistics must reproduce the calc_bunch_params values the tracker
     stored at element ends (proves the per-record datasets ARE sufficient).
  2. banked energy == the ledger's U_escaped column: same slices, independent
     bookkeeping paths (drain-time wavefront_params vs the zero-fill energy sum).
  3. analytic vs numerical propagation: each banked slice's rms size at the exit from
     the free-space ABCD map on its bank-time moment matrix must match the rms of the
     FFT-propagated slice in the pulse file (two independent routes).
  4. pulse pooling: the pooled whole-pulse sigma from per-slice params must match the
     directly computed moments of the concatenated pulse file.
  5. thread invariance: every dataset of stats.h5, escaped and pulse files identical
     at 1 vs 8 threads (dataset-level: HDF5 headers embed creation times).
  6. refusal: a dump_beam_at entry matching no element is refused by name.

Run by the benchmark harness; exits nonzero on failure.
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

import h5py
import numpy as np

from nml import to_groups

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from bunch_params_from_stats import bunch_params_at, pool_wavefront, M_ELECTRON  # noqa: E402
from read_stats import read_stats, same_data  # noqa: E402
from beamio import read_slices  # noqa: E402

import validate_bmad_stats  # noqa: E402  The standard's own conformance checker.

FAILED = False

NML = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "{lat}"
  out_root = "{root}"
  lambda0 = 1e-10
  beam_init%n_particle = 512
  beam_init%bunch_charge = 8.0e-15
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -4e-10
  beam_init%grid(3)%x_max = 4e-10
  beam_init%sig_pz = 8.8045e-5
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  nbins = 8
  seed_power = 1e4
  seed_waist_size = 30e-6
  grid_n_pts = 63
  grid_half_width = 2e-4
  window_length = 8e-10
  window_sample = 1
  shotnoise = T
  ran_seed = 999
  keep_escaped_field = T
  dump_beam_at = "UND"
  dump_field_at = "UND"
{extra}&end
"""

WRAP = """call, file = aramis_1seg.bmad
fel_unaveraged = 1
wiggler::*[FEL_TRACKING] = fel_unaveraged
"""

# Zero-length wake elements, both polarities, for check 7. A zero-length element carrying
# an sr wake is a standard Bmad idiom and legitimate: a wake assigned over an element
# range or a class lands on zero-length members as a matter of course. What differs is
# whether the kick can act. scale_with_length = F means it does, and the element must be
# tracked, take a record, and repeat coords/s at the plane the element before it ended
# on. T means the kick is identically zero (wake_mod scales by l$), so the element is
# skipped and the file must come out as though the wake were not there at all. WKT_NONE
# replaces the T pipes with plain ones to make that comparison.
ZL_WAKE = """no_digested
parameter[geometry] = open
parameter[particle] = electron
parameter[e_tot] = 11357.82 * m_electron
beginning[beta_a] = 15
beginning[beta_b] = 15
{wkt}
WKF: pipe, l = 0, sr_wake = {{amp_scale = 1, scale_with_length = F,
  longitudinal = {{1e12, 0, 0, 0.25, none}}}}
UND: wiggler, l = 0.45, l_period = 0.015, field_calc = helical_model, &
      b_max = 0.84853 * (twopi / 0.015) * m_electron / c_light, &
      tracking_method = custom, ds_step = 0.015
SEG: line = (WKT, UND, WKF, WKT)
use, SEG
"""

WKT_WAKE = """WKT: pipe, l = 0, sr_wake = {amp_scale = 1, scale_with_length = T,
  longitudinal = {1e12, 0, 0, 0.25, none}}"""

WKT_NONE = """WKT: pipe, l = 0"""


def h5_identical(fa, fb):
    """Whether two HDF5 files hold the same data, meta/ excluded.

    meta/ IS EXCLUDED, deliberately. Provenance moved from attributes to datasets for
    HDF5's 64 kB attribute cap (manual sec:meta), and input_echo carries out_root while
    timestamp carries the clock, so any two runs differ there by construction. Nothing in
    meta/ is physics. Before the move the exclusion existed only by the accident of being
    attributes, which is not something a reader of this check could have reasoned about.

    same_data counts NaN as equal to NaN, which the stats file needs: see read_stats.
    """
    with h5py.File(fa) as a, h5py.File(fb) as b:
        names_a, names_b = [], []
        keep = lambda n, o: isinstance(o, h5py.Dataset) and not n.startswith("meta/")
        a.visititems(lambda n, o: names_a.append(n) if keep(n, o) else None)
        b.visititems(lambda n, o: names_b.append(n) if keep(n, o) else None)
        if sorted(names_a) != sorted(names_b):
            return False
        return all(same_data(a[n][()], b[n][()]) for n in names_a)


def check(name, value, tol, note=""):
    global FAILED
    ok = value <= tol
    print(f"--- {name}: {value:.3e} (check {tol:.0e}) {note} {'ok' if ok else '** FAIL **'}")
    if not ok:
        FAILED = True


def run(exe, wd, name, text, threads="8"):
    (wd / (name + ".nml")).write_text(to_groups(text))
    r = subprocess.run([str(exe), name + ".nml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": threads, "PATH": "/usr/bin:/bin"})
    if r.returncode != 0:
        print(f"FAIL: {name} exited {r.returncode}:\n{r.stdout[-3000:]}\n{r.stderr[-1000:]}")
        sys.exit(1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", required=True)
    ap.add_argument("--latdir", required=True)
    ap.add_argument("--workdir", required=True)
    args = ap.parse_args()

    exe = pathlib.Path(args.exe).resolve()
    latdir = pathlib.Path(args.latdir).resolve()
    wd = pathlib.Path(args.workdir).resolve()
    wd.mkdir(parents=True, exist_ok=True)
    (wd / "aramis_1seg.bmad").write_bytes((latdir / "aramis_1seg.bmad").read_bytes())
    (wd / "dg_wrap.bmad").write_text(WRAP)

    run(exe, wd, "dg", NML.format(lat="dg_wrap.bmad", root="dg", extra=""))

    # 0. THE ACCEPTANCE TEST IS THE STANDARD'S OWN VALIDATOR, which knows no program:
    # a bmad-stats file conforms when it reports zero MUST failures, and everything it
    # checks is checked from the spec alone. The named checks below hold what it cannot
    # see: the physics identities and the declarations only this writer knows it owes.
    rep = validate_bmad_stats.validate(wd / "dg.stats.h5")
    for level, where, msg in rep.rows:
        if level == "MUST":
            print(f"    {level} {where}: {msg}")
    check("bmad-stats: the generic validator reports zero MUST failures", rep.n_must, 0.5,
          note=f"[{sum(1 for r in rep.rows if r[0] != 'INFO')} findings in "
               f"{len(rep.rows)} rows]")

    with read_stats(wd / "dg.stats.h5") as st:
        nslice = int(st.run["n_slice"])

        # 0a. THE FILE DESCRIBES ITSELF: a walk of every dataset, checking that it says
        # what it is and what its dimensions run over. This is what makes a generic
        # reader possible, so it is checked generically rather than name by name.
        bad = []
        for name, unit, label, axes, descrip in st.units_table():
            if unit == '' or descrip == '' or label == '':
                bad.append(f"{name}: no unit, label or description")
        for name, kind in st.kind.items():
            if kind == '' or st.group_description[name] == '':
                bad.append(f"{name}/: group carries no kind or description")
        check("stats: every dataset carries unit, label, description and axes", len(bad), 0.5,
              note=f"[{len(st.units_table())} datasets, {len(st.kind)} groups" +
                   (f"; first problem {bad[0]}]" if bad else "]"))

        # 0a2. THE ACCEPTANCE TEST IS A GENERIC LOAD: label every dimension of every
        # dataset from @axes alone. It fails if a name does not resolve to a coordinate,
        # if a dimension has no name, if a length disagrees with its coordinate, or if
        # one dataset names an axis twice, which is what would force a reader to guess
        # from a length or to dedupe by rule. Nothing here knows a dataset name.
        bad = []
        for name, unit, label, axes, descrip in st.units_table():
            shape = st[name].shape
            if len(axes) != len(shape):
                bad.append(f"{name}: {len(shape)} dimensions, {len(axes)} named")
                continue
            if len(set(axes)) != len(axes):
                bad.append(f"{name}: names an axis twice ({','.join(axes)})")
            for k, ax in enumerate(axes):
                if ax not in st.axis_names:
                    bad.append(f"{name}: @axes names {ax}, which is not an axis")
                elif len(st.coord(ax)) != shape[k]:
                    bad.append(f"{name}: dimension {k} is {shape[k]} long, {ax} is "
                               f"{len(st.coord(ax))}")
        check("stats: every dimension of every dataset labels itself", len(bad), 0.5,
              note=f"[{len(st.axis_names)} axes: {', '.join(st.axis_names)}" +
                   (f"; first problem {bad[0]}]" if bad else "]"))

        # 0a3. THE RECORD NUMBER IS THE AXIS, and s is a variable on it. s must be
        # non-decreasing, and where it repeats the two records must sit in DIFFERENT
        # elements, which is what a zero-length element applying a wake kick does: it
        # records at the plane the element before it ended on. Two records at one plane
        # INSIDE one element would be a defect in the walk that this axis choice hides.
        # The zero-length wake lattice below is what makes this fire on real duplicates.
        ds = np.diff(st.s)
        dup = np.flatnonzero(ds == 0)
        unexplained = [int(i) for i in dup if st.ix_ele[i] == st.ix_ele[i + 1]]
        check("stats: s is non-decreasing along the record axis", int(np.sum(ds < 0)), 0.5)
        check("stats: every repeated s straddles an element boundary", len(unexplained), 0.5,
              note=f"[{len(dup)} repeats in {len(st.s)} records]")

        # 0a4. PROVENANCE IS DATA, NOT AN ATTRIBUTE. HDF5 caps one attribute at 64 kB
        # (measured: 65495 bytes is the largest that writes), and the echoed namelist is
        # already 12 kB with a real lattice text at 37 kB, so meta/ as attributes was at
        # half the cap with a warning for a failure path: the provenance would have gone
        # missing silently. Every text in meta/ is a dataset now, and nothing anywhere in
        # the file leans on the cap.
        meta = st.meta
        big = []

        def note_attrs(name, obj):
            for key, val in obj.attrs.items():
                if isinstance(val, bytes) and len(val) > 60000:
                    big.append(f"{name}@{key} is {len(val)} bytes")

        with h5py.File(wd / "dg.stats.h5") as h5:
            h5.visititems(note_attrs)
            note_attrs("/", h5)
            n_meta_attrs = len([k for k in h5["meta"].attrs if k not in ("kind", "description")])
        check("stats: no attribute is near HDF5's 64 kB cap", len(big), 0.5,
              note=f"[{'; '.join(big) if big else 'largest well under'}]")
        check("stats: meta/ holds its texts as datasets", n_meta_attrs, 0.5,
              note=f"[{len(meta)} datasets, input_echo {len(meta['input_echo'])} bytes]")

        # 0a5. AND IT DOES NOT OVERSTATE ITSELF. dg tracks a WRAPPER lattice, so the
        # parser opened two files and lattice_source holds only the outer one. Recording
        # the count is what keeps that honest, and this run is the case that proves it:
        # a reader of lattice_source alone would think it had the lattice.
        check("stats: n_lattice_files reports the wrapper's second file",
              abs(int(meta["n_lattice_files"]) - 2), 0.5,
              note=f"[{meta['n_lattice_files']} files, source "
                   f"{len(meta['lattice_source'])} bytes]")

        # 0a6. A STATS FILE TRAVELS. Nothing identifies a person or a machine unless the
        # run asked for it, and a path is a basename.
        leaks = [k for k in ("user", "cwd") if k in meta]
        leaks += ["lattice_file has a directory"] if "/" in meta["lattice_file"] else []
        check("stats: the default file carries no machine-local values", len(leaks), 0.5,
              note=f"[{', '.join(leaks) if leaks else 'clean'}]")

        # 0a7. THE SLICE COORDINATES ARE EXACT, AND ONE OF THEM IS NOT. The grid is
        # uniform in TIME: the migration invariant carries each particle's own beta, and
        # z = -beta*c*(t-t_ref), so at a grid point the beta CANCELS and the arrival-time
        # separation is ct_slice/c with no beta anywhere. t_slice must therefore be
        # -ct_slice/c EXACTLY, which it was not until 2.3 (it carried a spurious beta0,
        # 3.9e-9). z_slice is the one that needs a reference, so it is only equal to
        # beta0*ct_slice to the rounding of two ways of forming beta0.
        c_light = 299792458.0
        p0c = float(st.run["p0c"])
        beta0 = p0c / np.sqrt(p0c**2 + M_ELECTRON**2)
        ct = st.ct_slice
        check("stats: t_slice is -ct_slice/c exactly (abs)",
              float(np.max(np.abs(st.t_slice + ct / c_light))), 1e-30)
        rel = np.abs(st.z_slice[1:] / (beta0 * ct[1:]) - 1)
        check("stats: z_slice is beta0*ct_slice", float(np.max(rel)), 1e-15,
              note=f"[beta0 = 1 - {1 - beta0:.3e}]")
        check("stats: the slice axis is the index", abs(len(st.slice) - nslice), 0.5,
              note=f"[{nslice} slices, positions ct_slice, t_slice, z_slice on it]")

        # 0a8. THE ENVELOPE DATA ARE THE PARTICLES'. rel_max/rel_min are order
        # statistics relative to the stored centroid, and the dump at the UND end holds
        # the SAME particles the record saw, so the position entries must match to the
        # bit: the file stores x verbatim, and IEEE subtraction of the same centroid is
        # deterministic. The momentum entries cross the dump's unit round trip
        # (px*p0c on write, /p0c on read), so they get an ulp-scale tolerance.
        irec = int(np.flatnonzero(st.at_end)[0])
        rel_hi = st["beam/slice/rel_max"][irec]
        rel_lo = st["beam/slice/rel_min"][irec]
        cen_r = st["beam/slice/centroid"][irec]
        dump_sl = read_slices(wd / "dg-at1-UND.beam.h5", wavelength=1e-10)
        p0_mc_l = float(st.run["p0c"]) / M_ELECTRON
        worst_pos, worst_mom = 0.0, 0.0
        for isl, sd in enumerate(dump_sl):
            if sd["n"] == 0:
                continue
            for j, arr in ((0, sd["x"]), (2, sd["y"])):
                worst_pos = max(worst_pos,
                                abs(np.max(arr - cen_r[isl, j]) - rel_hi[isl, j]),
                                abs(np.min(arr - cen_r[isl, j]) - rel_lo[isl, j]))
            for j, arr in ((1, sd["px"] / p0_mc_l), (3, sd["py"] / p0_mc_l),
                           (5, sd["pz"])):
                sc = max(abs(rel_hi[isl, j]), abs(rel_lo[isl, j]), 1e-30)
                worst_mom = max(worst_mom,
                                abs(np.max(arr - cen_r[isl, j]) - rel_hi[isl, j]) / sc,
                                abs(np.min(arr - cen_r[isl, j]) - rel_lo[isl, j]) / sc)
        check("envelope: position extremes vs the dump's particles, exactly (abs)",
              worst_pos, 1e-30, note=f"[{len(dump_sl)} slices, stored centroid]")
        check("envelope: momentum extremes across the dump's unit round trip",
              worst_mom, 1e-12)

        # And the time entry is the z entry through the exact map t - <t> =
        # -(z - <z>)/(beta0 c): the largest time offset is the smallest z one.
        beta0_l = p0_mc_l / np.sqrt(p0_mc_l**2 + 1)
        c_l = 299792458.0
        rmx = st["beam/slice/rel_max"]
        rmn = st["beam/slice/rel_min"]
        fin = np.isfinite(rmx[..., 6])
        worst_t = float(np.max(np.abs(rmx[..., 6][fin] +
                                      rmn[..., 4][fin] / (beta0_l * c_l))))
        worst_t = max(worst_t, float(np.max(np.abs(rmn[..., 6][fin] +
                                                   rmx[..., 4][fin] / (beta0_l * c_l)))))
        check("envelope: the t entry is the z entry through -dz/(beta0 c), exactly (abs)",
              worst_t, 1e-30)

        # 0b. THE JOIN. at_element_end selects the element-end rows, and every record
        # sits inside the element its ix_ele names. An off-by-one in either would break
        # every layout plot and no physics check would see it.
        n_end = int(st.at_end.sum())
        check("stats: at_element_end selects exactly the element ends",
              abs(n_end - int(st.run["n_element_end"])), 0.5,
              note=f"[{n_end} of {len(st.s)} records]")
        check("stats: element-end arrays are aligned with the mask",
              abs(st["beam/bunch/centroid"].shape[0] - n_end), 0.5)
        s1, s2 = st["lattice/s_start"], st["lattice/s_end"]
        ix = st.ix_ele
        outside = int(np.sum((st.s < s1[ix] - 1e-9) | (st.s > s2[ix] + 1e-9)))
        check("stats: every record sits inside the element ix_ele names", outside, 0.5,
              note=f"[{len(st.s)} records against {len(s1)} lattice rows]")

        # 1. Reconstruction at the element end (= the last record). The entry is an index
        # into the group's own axis, looked up rather than assumed: that lookup is the
        # only thing standing between a reader and a silent transposition. The projected
        # planes and the normal modes are separate groups over separate axes, because an
        # eigen-emittance and a projected emittance are different quantities.
        plane = {name: i for i, name in enumerate(st.coord("plane"))}
        mode = {name: i for i, name in enumerate(st.coord("mode"))}
        worst = 0.0
        for isl in range(nslice):
            bp = bunch_params_at(st, -1, isl)
            for m in ("x", "y", "z"):
                for pn in ("beta", "alpha", "emit", "norm_emit", "sigma", "sigma_p", "eta", "etap"):
                    stored = st[f"beam/slice_twiss/twiss/{pn}"][0, isl, plane[m]]
                    worst = max(worst, abs(bp[m][pn] - stored) / max(abs(stored), 1e-30))
            # Bmad's mat_eigen labels modes by eigenvector structure rather than by
            # magnitude, which coords/mode says, so this compares them as a sorted set.
            modes = sorted(st["beam/slice_twiss/modes/emit"][0, isl, mode[m]]
                           for m in ("a", "b", "c"))
            mine = sorted(bp[m]["emit"] for m in ("a", "b", "c"))
            for s0, m0 in zip(modes, mine):
                worst = max(worst, abs(m0 - s0) / max(abs(s0), 1e-30))
        # Measured 4.6e-8 (numpy eig vs mat_eigen). The projected planes agree exactly.
        check("stats: bunch_params reconstruction vs stored calc_bunch_params", worst, 1e-6)

        # 1b. The whole-window row is ASSEMBLED from the per-slice moments by the
        # pooled-covariance identity rather than summed over particles. Re-implement the
        # identity here, independently, including the local-to-global z map: a shift of
        # beta*(is-1)*spacing plus the pz shear that map carries because beta is the
        # particle's own. Measured against the particle sum when this landed: 4.0e-12 on
        # this config, 5.0e-11 on the 96-slice SASE example.
        spacing = float(st.run["slice_spacing"])
        p0_mc = float(st.run["p0c"]) / M_ELECTRON
        worst_pool = 0.0
        worst_nobg = 0.0
        for ie in range(st["beam/bunch/sigma"].shape[0]):
            cen = st["beam/slice_twiss/centroid"][ie]
            sig = st["beam/slice_twiss/sigma"][ie]
            wsl = st["beam/slice_twiss/charge_live"][ie]
            wtot = wsl.sum()
            if wtot <= 0:
                continue
            pmc = p0_mc * (1.0 + cen[:, 5])
            beta = pmc / np.sqrt(pmc**2 + 1)
            ell = np.arange(len(wsl)) * spacing
            shear = ell * p0_mc / np.sqrt(pmc**2 + 1)**3
            cg = cen.copy()
            cg[:, 4] = cg[:, 4] + beta * ell
            m = (wsl[:, None] * cg).sum(0) / wtot
            pooled = np.zeros((6, 6))
            nobg = np.zeros((6, 6))
            for isl in range(len(wsl)):
                s6 = sig[isl].copy()
                k = shear[isl]
                s0 = sig[isl]
                s6[4, 4] = s0[4, 4] + 2 * k * s0[4, 5] + k * k * s0[5, 5]
                for j in range(6):
                    if j == 4:
                        continue
                    s6[4, j] = s0[4, j] + k * s0[5, j]
                    s6[j, 4] = s6[4, j]
                d = cg[isl] - m
                pooled += wsl[isl] * (s6 + np.outer(d, d))
                nobg += wsl[isl] * s6
            pooled /= wtot
            nobg /= wtot
            st_c = st["beam/bunch/centroid"][ie]
            st_s = st["beam/bunch/sigma"][ie]
            den = np.maximum(np.abs(st_s), np.abs(pooled))
            worst_pool = max(worst_pool, float(np.max(np.where(den > 0, np.abs(st_s - pooled) / np.where(den > 0, den, 1), 0.0))))
            dc = np.maximum(np.abs(st_c), np.abs(m))
            worst_pool = max(worst_pool, float(np.max(np.where(dc > 0, np.abs(st_c - m) / np.where(dc > 0, dc, 1), 0.0))))
            den = np.maximum(np.abs(st_s), np.abs(nobg))
            worst_nobg = max(worst_nobg, float(np.max(np.where(den > 0, np.abs(st_s - nobg) / np.where(den > 0, den, 1), 0.0))))
        check("stats: whole-window row vs the pooled-covariance identity", worst_pool, 1e-9)
        # The check has teeth: dropping the between-group term (m_s - m)(m_s - m)^T --
        # the term that makes the identity exact -- must be caught, not absorbed. It
        # collapses sigma(5,5) from the window length squared to a slice length squared.
        check("stats: the between-group term is load-bearing (0 = confirmed)",
              0.0 if worst_nobg > 1e-3 else 1.0, 0.5,
              note=f"[dropping it moves the row by {worst_nobg:.2e}]")

        # 1c. The moments the element-end twiss rests on are the RECORD's, not a second
        # copy: beam/slice_twiss/centroid at an element end must be beam/slice/centroid
        # at the record the mask selects. Exactly, since one is written from the other.
        worst_ef = float(np.max(np.abs(st["beam/slice/centroid"][st.at_end] -
                                      st["beam/slice_twiss/centroid"])))
        worst_ef = max(worst_ef, float(np.max(np.abs(st["beam/slice/sigma"][st.at_end] -
                                                    st["beam/slice_twiss/sigma"]))))
        check("stats: the element-end moments ARE the record's (abs)", worst_ef, 1e-30)

        # And the field's theta moments are computed where the file says they are:
        # (nz, ns) per slice, at every element end that HAS field, plus the initial
        # record (the walk takes that one with angles so the starting state is
        # complete). An empty slice has no moments to compute and says so.
        valid = st["field/x/angle_moments_valid"].astype(bool)
        has_field = st["field/x/power"] > 0
        rows = np.zeros(len(st.s), bool)
        rows[0] = True
        rows |= st.at_end
        want = has_field & rows[:, None]
        check("stats: angle moments valid exactly where they were computed",
              int(np.sum(valid != want)), 0.5,
              note=f"[{int(valid.sum())} of {valid.size} slice-records]")

        # NaN, not a zero that reads as an answer, everywhere else. This is what lets a
        # consumer tell an empty slice's missing moments from a real measurement.
        cen = st["field/x/centroid"]
        check("stats: theta moments are NaN where they were not computed",
              0.0 if np.all(np.isnan(cen[..., 1][~valid])) else 1.0, 0.5)
        check("stats: and finite where they were",
              0.0 if np.all(np.isfinite(cen[..., 1][valid])) else 1.0, 0.5)

        f_cen = st["field/x/centroid"][-1]
        f_sig = st["field/x/sigma"][-1]
        f_en = st["field/total/energy"][-1]

    # 2. Banked energy == ledger U_escaped.
    led = np.loadtxt(wd / "dg.ledger.txt")
    with h5py.File(wd / "dg-escaped.fld.h5") as h5:
        nb = int(h5["slicecount"][0])
        pms = np.array([h5[f"slice{k:06d}/wavefront_params"][:] for k in range(1, nb + 1)])
        zt = np.array([h5[f"slice{k:06d}/z_transmit"][0] for k in range(1, nb + 1)])
    e_banked = pms[:, 20].sum()
    check("bank: banked slice energies vs ledger U_escaped, rel",
          abs(e_banked / led[-1, 4] - 1), 1e-12)

    # 3. Analytic (moment-map) vs numerical (FFT) propagation to the exit plane.
    with h5py.File(wd / "dg-pulse.fld.h5") as h5:
        n_all = int(h5["slicecount"][0])
        nx = int(h5["gridpoints"][0])
        dx = float(h5["gridsize"][0])
        nlive = n_all - nb
        ax = (np.arange(nx) - 0.5 * (nx - 1)) * dx
        sxx_num, syy_num, w_num = [], [], []
        cx_num, cy_num = [], []
        for j in range(nlive + 1, n_all + 1):
            g = h5[f"slice{j:06d}"]
            fr = g["field-real"][:].reshape(nx, nx)   # x fastest: columns = x.
            fi = g["field-imag"][:].reshape(nx, nx)
            I = fr**2 + fi**2
            w = I.sum()
            w_num.append(w)
            if w == 0:
                sxx_num.append(0); syy_num.append(0); cx_num.append(0); cy_num.append(0)
                continue
            cx = (I * ax[None, :]).sum() / w
            cy = (I * ax[:, None]).sum() / w
            cx_num.append(cx); cy_num.append(cy)
            sxx_num.append((I * (ax[None, :] - cx) ** 2).sum() / w)
            syy_num.append((I * (ax[:, None] - cy) ** 2).sum() / w)
    # Pulse index nlive+j holds banked slice nb-j+1: reverse into transmission order.
    sxx_num = np.array(sxx_num)[::-1]
    syy_num = np.array(syy_num)[::-1]
    w_num = np.array(w_num)[::-1]

    z_end = led[-1, 0]  # Ledger rows end at the segment exit for this one-element line.
    dz = z_end - zt
    # Flattened 4x4 column-major: (i,j) -> 4*(j-1)+(i-1).
    s11, s21, s22 = pms[:, 4 + 0], pms[:, 4 + 1], pms[:, 4 + 5]
    s33, s43, s44 = pms[:, 4 + 10], pms[:, 4 + 11], pms[:, 4 + 15]
    sxx_ana = s11 + 2 * dz * s21 + dz**2 * s22
    syy_ana = s33 + 2 * dz * s43 + dz**2 * s44
    m = w_num > 1e-6 * w_num.max()
    rel = np.abs(np.sqrt(sxx_num[m] / sxx_ana[m]) - 1)
    rel = np.maximum(rel, np.abs(np.sqrt(syy_num[m] / syy_ana[m]) - 1))
    # Measured max 8.4e-3 (grid-edge wrap of the FFT route on divergent noise slices).
    check("bank: analytic vs FFT-propagated rms at exit, max rel", rel.max(), 2e-2)

    # 4. Pulse pooling: pooled per-slice params vs directly computed pulse moments.
    #    Banked side: params propagated analytically to the exit. Live side: the last
    #    stats record (angles valid there). Compare sigma_xx of the pooled result
    #    against pooling the pulse file's numerically computed per-slice moments.
    cen_b = np.stack([pms[:, 0] + dz * pms[:, 1], pms[:, 1],
                      pms[:, 2] + dz * pms[:, 3], pms[:, 3]], axis=1)
    sig_b = np.zeros((nb, 4, 4))
    sig_b[:, 0, 0] = sxx_ana;  sig_b[:, 2, 2] = syy_ana
    en_b = pms[:, 20]
    cen_l = f_cen;  sig_l = f_sig;  en_l = f_en
    _, pooled_ana, _ = pool_wavefront(np.vstack([cen_l, cen_b]),
                                      np.vstack([sig_l.reshape(-1, 16), sig_b.reshape(-1, 16)]),
                                      np.concatenate([en_l, en_b]))
    # Numerical pooling over the SAME banked slices from the pulse file:
    cen_n = np.stack([np.array(cx_num)[::-1], np.zeros(nb), np.array(cy_num)[::-1], np.zeros(nb)], axis=1)
    sig_n = np.zeros((nb, 4, 4))
    sig_n[:, 0, 0] = sxx_num;  sig_n[:, 2, 2] = syy_num
    en_n = w_num * (en_b.sum() / max(w_num.sum(), 1e-300))   # Same total, per-slice proportional.
    _, pooled_num, _ = pool_wavefront(np.vstack([cen_l, cen_n]),
                                      np.vstack([sig_l.reshape(-1, 16), sig_n.reshape(-1, 16)]),
                                      np.concatenate([en_l, en_n]))
    check("pulse: pooled sigma_xx analytic vs numerical routes, rel",
          abs(pooled_ana[0, 0] / pooled_num[0, 0] - 1), 2e-2)

    # 5. Thread invariance: every dataset of all three files identical at 1 vs 8
    #    threads. Dataset-level, not raw bytes: HDF5 object headers embed creation
    #    times, so two files with bit-identical DATA differ as byte streams -- the
    #    data is the invariance claim, the container metadata is not.
    run(exe, wd, "dg1", NML.format(lat="dg_wrap.bmad", root="dg1", extra=""), threads="1")

    same = all(h5_identical(wd / f"dg{s}", wd / f"dg1{s}")
               for s in (".stats.h5", "-escaped.fld.h5", "-pulse.fld.h5"))
    check("thread invariance: stats/escaped/pulse data identical 1 vs 8 (1 = yes)",
          0.0 if same else 1.0, 0.5)

    # 5a. ZERO-LENGTH WAKE ELEMENTS, BOTH POLARITIES. This is the lattice the repeated-s
    # check was written for, and until it existed that check had only ever seen zero
    # repeats, which is to say it was untested. WKF can kick, so it is tracked and its
    # record lands on the plane UND ended at, repeating coords/s. The two WKT pipes
    # cannot kick, so they are skipped and the file must be IDENTICAL to one whose
    # zero-length pipes carry no wake at all.
    (wd / "zlw.bmad").write_text(ZL_WAKE.format(wkt=WKT_WAKE))
    (wd / "zlwn.bmad").write_text(ZL_WAKE.format(wkt=WKT_NONE))
    run(exe, wd, "zlw", NML.format(lat="zlw.bmad", root="zlw", extra=""), threads="4")
    run(exe, wd, "zlwn", NML.format(lat="zlwn.bmad", root="zlwn", extra=""), threads="4")

    with read_stats(wd / "zlw.stats.h5") as st:
        ds = np.diff(st.s)
        dup = np.flatnonzero(ds == 0)
        unexplained = [int(i) for i in dup if st.ix_ele[i] == st.ix_ele[i + 1]]
        names = st.ele_name
        check("zero-length wake: coords/s repeats where one CAN kick", 1.0 / max(len(dup), 1) - 1.0,
              0.5, note=f"[{len(dup)} repeat(s), at " +
                        ", ".join(f"{names[i]}->{names[i+1]}" for i in dup) + "]")
        check("zero-length wake: and every repeat straddles an element boundary",
              len(unexplained), 0.5)
        check("zero-length wake: s is still non-decreasing", int(np.sum(ds < 0)), 0.5)

    same = h5_identical(wd / "zlw.stats.h5", wd / "zlwn.stats.h5")
    check("zero-length wake: one that CANNOT kick leaves the file unchanged (0 = yes)",
          0.0 if same else 1.0, 0.5)

    # 5b. The environment switch: opt-in, so the fields appear only when asked for.
    run(exe, wd, "dgenv", NML.format(lat="dg_wrap.bmad", root="dgenv",
        extra="  global%record_environment = T\n"), threads="2")
    with read_stats(wd / "dgenv.stats.h5") as st:
        got = sorted(k for k in st.meta if k in ("user", "cwd"))
    check("provenance: global%record_environment restores user and cwd",
          abs(len(got) - 2), 0.5, note=f"[{', '.join(got) if got else 'neither'}]")

    # 6. Refusal: dump list entry matching nothing, refused by name.
    (wd / "dg_bad.nml").write_text(to_groups(NML.format(lat="dg_wrap.bmad", root="dg_bad",
                                             extra='  dump_beam_at(2) = "NO_SUCH_ELEMENT"\n')))
    r = subprocess.run([str(exe), "dg_bad.nml"], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
    refused = r.returncode != 0 and "NO_SUCH_ELEMENT" in (r.stdout + r.stderr)
    check("refusal: unknown dump_beam_at element refused by name (1 = yes)",
          0.0 if refused else 1.0, 0.5)

    if FAILED:
        print("DIAGNOSTIC CHECKS: FAIL")
        sys.exit(1)
    print("diagnostic checks: PASS")


if __name__ == "__main__":
    sys.exit(main())
