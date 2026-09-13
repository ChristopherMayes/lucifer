#!/usr/bin/env python3
"""
Self-referenced checks for the collective terms (the Genesis-comparison
tiers live in compare_fel.py).

1. Energy bookkeeping. A dark, quiet, cold beam (no seed, no shot noise, negligible
   energy spread) tracked with all wake kernels on has exactly one energy channel: the
   applied per-slice eloss. At every consecutive record pair, the measured
   d<gamma> must equal eloss(slice)*dz/m_electron to the diagnostic's print resolution,
   and sigma_gamma must stay constant (the kick is uniform within a slice).

2. Stale-wake structure. With migration on and a large energy spread, the current
   profile changes, and the wake convolution must follow: <out_root>.wake.txt must hold
   more than one z-stamped eloss block, and the blocks must differ. Removing the
   migration-stride recompute (the stale-wake mutation) leaves one block and fails
   loudly. The checker parses the driver's record of recomputes rather than
   reimplementing the convolution.

3. A stride that only drops. A particle dropped off the window's end changes the
   current profile as a move does, and the recompute once keyed on moves alone, so a
   stride that dropped charge and moved nothing left the wake computed from charge no
   longer there. The case is built rather than found: a cold beam confined to the head
   slice by editing a dump, detuned so it slips out of the window and never into
   another slice. Its migration record must show drops and no moves, and the wake
   record must hold one block per such stride beyond the hoist.

Usage: check_collective.py --exe <lucifer> --workdir <dir>
The workdir must hold aramis_1seg.bmad and aramis.bmad. Exit 0 only if all pass.
"""

from __future__ import annotations

import argparse
import math
import shutil
import pathlib
import re
import subprocess
import sys

import h5py
import numpy as np

from nml import to_groups

M_ELECTRON = 0.51099895069e6
EPS0 = 8.8541878128e-12
SPACING = 3e-10   # the probe deck's slice spacing, n_wavelength * lambda0

BASE = """! flat keys; routed into the three groups by nml.to_groups
  lat_file = "{lat}"
  out_root = "{root}"
  source_filter = F
  lambda0 = 1e-10
  beam_init%n_particle = 512
  beam_init%bunch_charge = 2.401661485427e-14
  beam_init%distribution_type(3) = "GRID"
  beam_init%grid(3)%x_min = -1.200000e-09
  beam_init%grid(3)%x_max = 1.200000e-09
  beam_init%sig_pz = {sig_pz}
  beam_init%a_norm_emit = 4e-7
  beam_init%b_norm_emit = 4e-7
  seed_power = 0
  grid_n_pts = 64
  grid_half_width = 2e-4
  beamlet_size = 8
  window_length = 2.4e-9
  n_wavelength = 3
  ran_seed = 555
  wake_on = T
  wake_radius = 2.5e-3
  wake_conductivity = 5.813e7
  wake_relaxation = 8.1e-6
  wake_gap = 0.5e-3
  wake_lgap = 0.015
  wake_hrough = 100e-9
  wake_lrough = 100e-6
  write_diag = T
{extra}&end
"""


def run(exe, wd, name, text):
    (wd / name).write_text(to_groups(text))
    r = subprocess.run([str(exe), name], cwd=wd, capture_output=True, text=True,
                       env={"OMP_NUM_THREADS": "8", "PATH": "/usr/bin:/bin"})
    if r.returncode != 0:
        print(f"FAIL: {name} exited {r.returncode}:\n{r.stdout[-2000:]}")
        sys.exit(1)
    return r.stdout


def migration_events(fn):
    """[(s, moved, charge_dropped)] per migration event, the per-event rows only."""
    events = []
    for line in pathlib.Path(fn).read_text().splitlines():
        f = line.split()
        if not f or line.startswith("#") or f[0] in ("moved", "charge_dropped_total",
                                                       "worst_bunching_deviation"):
            continue
        events.append((float(f[0]), int(f[1]), float(f[2])))
    return events


def split_particles(src, dst, parts, which=None):
    """Copy a beam dump, replacing particles by colocated copies of divided charge.

    The copies sit at the same coordinates and their weights sum to the original's, so
    the beam is the same physical charge distribution described with more macroparticles.
    Every field the slice produces, and every kick it receives, must be unchanged.

    which = None splits every particle, which preserves an unweighted centroid as well as
    a weighted one and is the control. An index splits that particle alone, which moves an
    unweighted centroid and leaves the weighted one where it was.
    """
    shutil.copy(src, dst)
    with h5py.File(dst, "r+") as h5:
        g = h5["data/00001/particles/electron"]
        n = g["particlePatches/numParticles"][()]
        off = g["particlePatches/numParticlesOffset"][()]
        n_tot = int(n.sum())
        take = np.arange(n_tot)
        reps = np.ones(n_tot, dtype=int)
        reps[take if which is None else [which]] = parts
        idx = np.repeat(take, reps)
        names = []
        g.visititems(lambda nm, o: names.append(nm)
                     if isinstance(o, h5py.Dataset) and o.shape and o.shape[0] == n_tot
                     and not nm.startswith("particlePatches") else None)
        for nm in names:
            data, attrs = g[nm][...][idx], dict(g[nm].attrs)
            del g[nm]
            d = g.create_dataset(nm, data=data)
            for key, val in attrs.items():
                d.attrs[key] = val

        # The charge each copy carries. The loop above has already repeated a shaped
        # weight record; a constant one becomes shaped here, since the copies no longer
        # carry what the rest of the beam carries.
        w = g["weight"]
        if isinstance(w, h5py.Dataset) and w.shape:
            wv = w[...] / reps[idx]
            w[...] = wv
        else:
            wv = np.full(len(idx), float(np.ravel(w.attrs["value"])[0])) / reps[idx]
            attrs = {k: v for k, v in w.attrs.items() if k not in ("value", "shape")}
            del g["weight"]
            d = g.create_dataset("weight", data=wv)
            for k, v in attrs.items():
                d.attrs[k] = v

        for nm, o in list(g.items()):
            if "shape" in o.attrs and int(np.ravel(o.attrs["shape"])[0]) == n_tot:
                o.attrs["shape"] = np.array([len(idx)], dtype=o.attrs["shape"].dtype)
        g.attrs["numParticles"] = np.array([len(idx)], dtype=g.attrs["numParticles"].dtype)
        new_n = np.array([int(reps[off[k]:off[k] + n[k]].sum()) for k in range(len(n))])
        g["particlePatches/numParticles"][...] = new_n
        g["particlePatches/numParticlesOffset"][...] = np.concatenate(
            [[0], np.cumsum(new_n)[:-1]]).astype(off.dtype)
    return len(idx)


def collapse_slice(src, dst, islice, keep, zero_charge=False):
    """Copy a beam dump, leaving `keep` particles of one slice at one transverse point.

    The kept particles take the slice's centroid in x, y, px and py, so they stay at one
    point through the element rather than re-spreading from their own transverse momenta,
    and they keep their own z, so they sit at different ponderomotive phases. keep = 1 is
    the state migration leaves behind. The charge is preserved, or zeroed on request.
    """
    shutil.copy(src, dst)
    with h5py.File(dst, "r+") as h5:
        g = h5["data/00001/particles/electron"]
        n = g["particlePatches/numParticles"][()]
        off = g["particlePatches/numParticlesOffset"][()]
        a, b = int(off[islice]), int(off[islice] + n[islice])
        n_tot = int(n.sum())
        idx = np.concatenate([np.arange(0, a + keep), np.arange(b, n_tot)])
        names = []
        g.visititems(lambda nm, o: names.append(nm)
                     if isinstance(o, h5py.Dataset) and o.shape and o.shape[0] == n_tot
                     and not nm.startswith("particlePatches") else None)
        w = g["weight"]
        wv = (w[...] if (isinstance(w, h5py.Dataset) and w.shape)
              else np.full(n_tot, float(np.ravel(w.attrs["value"])[0])))
        w_slice = wv[a:b].sum()
        mid = {nm: float(g[nm][a:b].mean()) for nm in
               ("position/x", "position/y", "momentum/x", "momentum/y")}
        for nm in names:
            data, attrs = g[nm][...][idx], dict(g[nm].attrs)
            del g[nm]
            d = g.create_dataset(nm, data=data)
            for k, v in attrs.items():
                d.attrs[k] = v
        for nm, v in mid.items():
            d = g[nm]
            arr = d[...]
            arr[a:a + keep] = v            # one transverse point, one transverse momentum
            d[...] = arr
        w = g["weight"]
        if not (isinstance(w, h5py.Dataset) and w.shape):
            attrs = {k: v for k, v in w.attrs.items() if k not in ("value", "shape")}
            del g["weight"]
            d = g.create_dataset("weight", data=wv[idx])
            for k, v in attrs.items():
                d.attrs[k] = v
        g["weight"][a:a + keep] = 0.0 if zero_charge else w_slice / keep
        for nm, o in list(g.items()):
            if "shape" in o.attrs and int(np.ravel(o.attrs["shape"])[0]) == n_tot:
                o.attrs["shape"] = np.array([len(idx)], dtype=o.attrs["shape"].dtype)
        g.attrs["numParticles"] = np.array([len(idx)], dtype=g.attrs["numParticles"].dtype)
        new_n = n.copy()
        new_n[islice] = keep
        g["particlePatches/numParticles"][...] = new_n
        g["particlePatches/numParticlesOffset"][...] = np.concatenate(
            [[0], np.cumsum(new_n)[:-1]]).astype(off.dtype)


def narrow_slice(src, dst, islice, factor=None, singleton=False):
    """Copy a beam dump, narrowing one slice: its spread scaled, or one particle left.

    Scaling shrinks the transverse spread about the slice's own centroid and keeps every
    particle. singleton keeps one particle carrying the slice's whole charge, which is
    the state migration reaches by moving particles out one at a time and is a slice of
    exactly zero width for as long as it lasts. Either way the charge and its place are
    what they were, so the long-range kernel sees the same source with a smaller area.
    """
    shutil.copy(src, dst)
    with h5py.File(dst, "r+") as h5:
        g = h5["data/00001/particles/electron"]
        n = g["particlePatches/numParticles"][()]
        off = g["particlePatches/numParticlesOffset"][()]
        a, b = int(off[islice]), int(off[islice] + n[islice])
        if not singleton:
            for nm in ("position/x", "position/y"):
                d = g[nm]
                v = d[...]
                mid = v[a:b].mean()
                v[a:b] = mid + (v[a:b] - mid) * factor
                d[...] = v
            return
        n_tot = int(n.sum())
        keep = np.concatenate([np.arange(0, a + 1), np.arange(b, n_tot)])
        names = []
        g.visititems(lambda nm, o: names.append(nm)
                     if isinstance(o, h5py.Dataset) and o.shape and o.shape[0] == n_tot
                     and not nm.startswith("particlePatches") else None)
        w = g["weight"]
        wv = (w[...] if (isinstance(w, h5py.Dataset) and w.shape)
              else np.full(n_tot, float(np.ravel(w.attrs["value"])[0])))
        w_slice = wv[a:b].sum()
        for nm in names:
            data, attrs = g[nm][...][keep], dict(g[nm].attrs)
            del g[nm]
            d = g.create_dataset(nm, data=data)
            for k, v in attrs.items():
                d.attrs[k] = v
        w = g["weight"]
        if not (isinstance(w, h5py.Dataset) and w.shape):
            attrs = {k: v for k, v in w.attrs.items() if k not in ("value", "shape")}
            del g["weight"]
            d = g.create_dataset("weight", data=wv[keep])
            for k, v in attrs.items():
                d.attrs[k] = v
        g["weight"][a] = w_slice          # the kept particle carries the slice's charge
        for nm, o in list(g.items()):
            if "shape" in o.attrs and int(np.ravel(o.attrs["shape"])[0]) == n_tot:
                o.attrs["shape"] = np.array([len(keep)], dtype=o.attrs["shape"].dtype)
        g.attrs["numParticles"] = np.array([len(keep)], dtype=g.attrs["numParticles"].dtype)
        new_n = n.copy()
        new_n[islice] = 1
        g["particlePatches/numParticles"][...] = new_n
        g["particlePatches/numParticlesOffset"][...] = np.concatenate(
            [[0], np.cumsum(new_n)[:-1]]).astype(off.dtype)


def keep_one_patch(src, dst, k):
    """Copy a beam dump keeping the particles of patch k only (-1 is the last, the head).

    A dump is one particlePatch per slice, so a beam confined to one slice is the dump
    with the other patches emptied: every per-particle record cut to the patch's range,
    the constant records' shape and the species' numParticles set to the count kept, and
    the patch table rewritten with one patch full and the offsets at zero.
    """
    shutil.copy(src, dst)
    with h5py.File(dst, "r+") as h5:
        g = h5["data/00001/particles/electron"]
        n = g["particlePatches/numParticles"][()]
        off = g["particlePatches/numParticlesOffset"][()]
        k = k % len(n)
        n_tot, a, b = int(n.sum()), int(off[k]), int(off[k] + n[k])
        names = []
        g.visititems(lambda nm, o: names.append(nm)
                     if isinstance(o, h5py.Dataset) and o.shape and o.shape[0] == n_tot
                     and not nm.startswith("particlePatches") else None)
        for nm in names:
            data, attrs = g[nm][a:b], dict(g[nm].attrs)
            del g[nm]
            d = g.create_dataset(nm, data=data)
            for key, val in attrs.items():
                d.attrs[key] = val
        def cut_shape(nm, o):
            if "shape" in o.attrs and int(np.ravel(o.attrs["shape"])[0]) == n_tot:
                o.attrs["shape"] = np.array([b - a], dtype=o.attrs["shape"].dtype)
        g.visititems(cut_shape)
        g.attrs["numParticles"] = np.array([b - a], dtype=g.attrs["numParticles"].dtype)
        kept = np.zeros_like(n)
        kept[k] = b - a
        g["particlePatches/numParticles"][...] = kept
        g["particlePatches/numParticlesOffset"][...] = np.zeros_like(off)
    return b - a


def wake_blocks(fn):
    """[(z, eloss array)] per z-stamped block."""
    blocks, z, vals = [], None, []
    for line in pathlib.Path(fn).read_text().splitlines():
        m = re.match(r"# z =\s+(\S+)", line)
        if m:
            if z is not None:
                blocks.append((z, np.array(vals)))
            z, vals = float(m.group(1)), []
        elif line.strip() and not line.startswith("#"):
            vals.append(float(line.split()[1]))
    if z is not None:
        blocks.append((z, np.array(vals)))
    return blocks


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    p.add_argument("--exe", required=True)
    p.add_argument("--workdir", required=True)
    a = p.parse_args()
    wd = pathlib.Path(a.workdir)
    exe = pathlib.Path(a.exe).resolve()
    ok = True

    # 1. Energy bookkeeping on a cold dark beam, one segment.
    run(exe, wd, "colle.nml", BASE.format(lat="aramis_1seg.bmad", root="colle",
                                          sig_pz="8.804506566858e-8", extra=""))
    eloss = wake_blocks(wd / "colle.wake.txt")[0][1]
    d = np.loadtxt(wd / "colle.diag.txt")
    ns = int(d[:, 1].max())
    d = d.reshape(-1, ns, d.shape[1])
    # The diag's energy columns are Bmad-convention eV, so the bookkeeping identity
    # is d<E> = eloss*dz directly (eloss is eV/m).
    z, me, se = d[:, 0, 0], d[:, :, 6], d[:, :, 7]
    dz = np.diff(z)
    de_meas = np.diff(me, axis=0)
    de_exp = eloss[None, :] * dz[:, None]
    err = np.abs(de_meas - de_exp).max()
    scale = np.abs(de_exp).max()
    e_ok = err < 1e-6 * scale + 1e-7 * M_ELECTRON
    ok = ok and e_ok
    print(f"--- collective energy bookkeeping: max |d<E> - eloss*dz| = {err:.2e} eV"
          f" against kicks of {scale:.2e} eV  {'ok' if e_ok else 'FAIL'}")
    s_dev = np.abs(np.diff(se, axis=0)).max()
    s_ok = s_dev < 1e-8 * M_ELECTRON
    ok = ok and s_ok
    print(f"--- collective sigma_energy invariance under uniform kicks: {s_dev:.2e} eV  "
          f"{'ok' if s_ok else 'FAIL'}")

    # 2. Stale-wake structure: heavy migration must force recomputes that change eloss.
    run(exe, wd, "collm.nml", BASE.format(lat="aramis.bmad", root="collm",
        sig_pz="5.282703940115e-3", extra="  migrate = T\n"))
    blocks = wake_blocks(wd / "collm.wake.txt")
    changed = len(blocks) > 1 and any(
        not np.array_equal(blocks[0][1], b[1]) for _, b in [(None, blk) for _, blk in blocks[1:]])
    m_ok = len(blocks) > 1 and changed
    ok = ok and m_ok
    print(f"--- stale-wake structure: {len(blocks)} eloss blocks under heavy migration, "
          f"changed = {changed}  {'ok' if m_ok else 'FAIL'}")

    # 3. A stride that only drops. The source run writes the initial state at the detuned
    # wavelength, the dump is cut to its head slice, and the cut beam is loaded back with
    # the field the source wrote. Two segments: the first stride moves nothing, the
    # second drops most of the slice off the head and moves nothing, and the rest of the
    # beam is still there for the final dump. A cold beam at a low charge, so the slip is
    # the detuning's alone and no gain spreads it into both directions.
    drop_extra = ("  write_initial = T\n  load_only = T\n")
    src = BASE.format(lat="aramis.bmad", root="colld_src", sig_pz="1e-8", extra=drop_extra)
    src = src.replace("lambda0 = 1e-10", "lambda0 = 1.006e-10").replace(
        "bunch_charge = 2.401661485427e-14", "bunch_charge = 2.401661485427e-16")
    run(exe, wd, "colld_src.nml", src)
    n_head = keep_one_patch(wd / "colld_src-initial.beam.h5", wd / "colld_head.beam.h5", -1)
    text = BASE.format(lat="aramis.bmad", root="colld", sig_pz="1e-8",
                       extra='  migrate = T\n  beam_file = "colld_head.beam.h5"\n'
                             '  field_file = "colld_src-initial.wf.h5"\n  load_mode = "keep"\n'
                             '  global%track_end = "UND##2"\n')
    text = text.replace("lambda0 = 1e-10", "lambda0 = 1.006e-10").replace(
        "bunch_charge = 2.401661485427e-14", "bunch_charge = 2.401661485427e-16")
    run(exe, wd, "colld.nml", text)
    events = migration_events(wd / "colld.migration.txt")
    drop_only = [e for e in events if e[1] == 0 and e[2] > 0]
    with_moves = [e for e in events if e[1] > 0]
    blocks_d = wake_blocks(wd / "colld.wake.txt")
    refreshed = len(blocks_d) == 1 + len(events) and all(
        not np.array_equal(blocks_d[0][1], b) for _, b in blocks_d[1:])
    d_ok = len(drop_only) >= 1 and not with_moves and refreshed
    ok = ok and d_ok
    print(f"--- drop-only stride: {n_head} particles in the head slice, {len(drop_only)} stride(s) "
          f"dropped and moved nothing, {len(with_moves)} moved; {len(blocks_d)} eloss blocks "
          f"for {len(events)} stride(s), refreshed = {refreshed}  {'ok' if d_ok else 'FAIL'}")

    # 4. The charge representation. Replacing a particle by colocated copies whose
    # charges sum to the original's is the same physical beam described with more
    # macroparticles, so every field it makes and every kick it takes must be unchanged.
    # The short-range solve centers its radial bins and its azimuthal basis on the slice
    # centroid, and that centroid was an unweighted mean where every source term is
    # charge weighted: splitting one particle moved the origin and changed the field at
    # every physical point. Splitting every particle equally leaves an unweighted mean
    # where it was, which is why the harness's own beamlet loader could not see it, so
    # the check splits one particle and keeps the all-particle case as the control.
    sc_lat = (wd / "collsc.bmad")
    sc_lat.write_text("call, file = aramis_1seg.bmad\nwiggler::*[SPACE_CHARGE_METHOD] = slice\n")
    sc_on = ("  bmad_com%csr_and_space_charge_on = T\n  space_charge%nz = 2\n"
             "  space_charge%nphi = 1\n  space_charge%rmax = 250e-6\n")
    # Seeded, since a dark run's field is the numerical floor and comparing two floors
    # says nothing about the solve that feeds them.
    seeded = BASE.replace("seed_power = 0", "seed_power = 1e4\n  seed_waist_size = 30e-6")
    src = seeded.format(lat="collsc.bmad", root="collsc_src", sig_pz="8.804506566858e-08",
                        extra="  write_initial = T\n  load_only = T\n")
    run(exe, wd, "collsc_src.nml", src)
    n_sel = split_particles(wd / "collsc_src-initial.beam.h5", wd / "collsc_sel.beam.h5",
                            16, which=0)
    n_all = split_particles(wd / "collsc_src-initial.beam.h5", wd / "collsc_all.beam.h5", 2)
    load = ('  beam_file = "{f}"\n  field_file = "collsc_src-initial.wf.h5"\n'
            '  load_mode = "keep"\n')
    for tag, extra in (("on", sc_on), ("off", "")):
        base_run = f"collsc_{tag}"
        for name, f in (("ref", "collsc_src-initial.beam.h5"), ("sel", "collsc_sel.beam.h5"),
                        ("all", "collsc_all.beam.h5")):
            run(exe, wd, f"{base_run}_{name}.nml",
                seeded.format(lat="collsc.bmad", root=f"{base_run}_{name}",
                              sig_pz="8.804506566858e-08", extra=load.format(f=f) + extra))
        ref = np.loadtxt(wd / f"{base_run}_ref.diag.txt")
        # The physical columns: power, on-axis intensity, bunching, the energy moments,
        # the beam sizes and the current. n_eff counts macroparticles and is meant to
        # move, and the bunching phase is an angle that wraps, so a difference across the
        # wrap says nothing that the bunching magnitude beside it does not say better.
        cols = [2, 3, 4, 6, 7, 8, 9, 10]
        # Each column against its own range over the run, so a quantity that passes
        # through zero is not compared against zero.
        scale = np.maximum(np.abs(ref).max(axis=0), 1e-30)
        worst = {}
        for name in ("sel", "all"):
            d = np.loadtxt(wd / f"{base_run}_{name}.diag.txt")
            worst[name] = float(np.max(np.abs(d[:, cols] - ref[:, cols]) / scale[cols]))
        s_ok = worst["sel"] < 1e-12 and worst["all"] < 1e-12
        ok = ok and s_ok
        print(f"--- charge representation, space charge {tag}: one particle split 16 ways "
              f"{worst['sel']:.2e}, every particle split {worst['all']:.2e} "
              f"(tol 1e-12)  {'ok' if s_ok else 'FAIL'}")

    # 5. The long-range kernel at a vanishing source width, against the closed form it
    # becomes there. A source slice contributes (1 - |d|/sqrt(d^2 + A)) / A with A its
    # transverse area and d the boosted separation, and that has a finite limit as A goes
    # to zero: 1 / 2 d^2, which carries the whole kernel to q / (4 pi eps0 d^2), the field
    # of a point charge. Written with the subtraction it loses every digit once A is small
    # against d^2 and divides by a zero A, and a guard put a square metre of area in its
    # place, so the field of a slice whose charge sits at one transverse point came out
    # orders too small. Migration reaches that state by moving particles out of a slice
    # one at a time until one is left, which is the deck here.
    #
    # The deck isolates the term: dark, so no gain, the chamber wake off and the radial
    # solve off, so a slice's mean energy moves only by the long-range term's work. The
    # source's own contribution is isolated by differencing against the same beam with
    # that slice's charge set to zero, which leaves every other source where it was.
    lr_only = ("  bmad_com%csr_and_space_charge_on = T\n  space_charge%nz = 0\n"
               "  space_charge%longrange = T\n  space_charge%rmax = 250e-6\n")
    dark = BASE.replace("wake_on = T", "wake_on = F")
    narrow_slice(wd / "collsc_src-initial.beam.h5", wd / "colllr_one.beam.h5", 0, singleton=True)
    narrow_slice(wd / "collsc_src-initial.beam.h5", wd / "colllr_nq.beam.h5", 0, singleton=True)
    with h5py.File(wd / "colllr_one.beam.h5") as h5:
        g = h5["data/00001/particles/electron"]
        q_one = float(g["weight"][int(g["particlePatches/numParticlesOffset"][()][0])])
    with h5py.File(wd / "colllr_nq.beam.h5", "r+") as h5:
        g = h5["data/00001/particles/electron"]
        g["weight"][int(g["particlePatches/numParticlesOffset"][()][0])] = 0.0
    got = {}
    for tag in ("one", "nq"):
        run(exe, wd, f"colllr_{tag}.nml",
            dark.format(lat="collsc.bmad", root=f"colllr_{tag}", sig_pz="8.804506566858e-08",
                        extra=load.format(f=f"colllr_{tag}.beam.h5") + lr_only))
        d = np.loadtxt(wd / f"colllr_{tag}.diag.txt")
        ns = int(d[:, 1].max())
        d = d.reshape(-1, ns, d.shape[1])
        got[tag] = (d[-1, :, 6] - d[0, :, 6], float(d[-1, 0, 0] - d[0, 0, 0]))
    dE, length = got["one"][0] - got["nq"][0], got["one"][1]
    # aw of the probe's wiggler, and the boost the kernel carries.
    gamma_z = 11357.82 / math.sqrt(1 + 0.84853**2)
    worst = 0.0
    for i in range(1, 4):                      # the three nearest, where the term dominates
        d_sep = i * SPACING * gamma_z
        want = q_one / (4 * math.pi * EPS0 * d_sep**2) * length
        worst = max(worst, abs(abs(dE[i]) / want - 1))
    p_ok = worst < 5e-2
    ok = ok and p_ok
    print(f"--- long-range kernel at zero width: a slice migration left one particle in, "
          f"its field on the three nearest slices against q / 4 pi eps0 d^2, worst "
          f"{worst:.2e} (tol 5e-2)  {'ok' if p_ok else 'FAIL'}")

    # 6. What the short-range solve does with the states a run can reach. Its grid needs a
    # radial scale, and a slice whose charge sits at one transverse point offers none: the
    # cell volumes go as the cell width squared, so the field has no limit there. That is
    # not a zero. Particles at one transverse point and different ponderomotive phases
    # push each other, and only a lone particle feels nothing, its own harmonic sum
    # cancelling. So the solve refuses where it has no scale, and where it has one it must
    # give the lone particle nothing and the pair their mutual push.
    sr_only = ("  bmad_com%csr_and_space_charge_on = T\n  space_charge%nz = 2\n"
               "  space_charge%nphi = 1\n")
    dark = BASE.replace("wake_on = T", "wake_on = F")

    def sc_run(root, beam_file, extra, expect_fail=False):
        text = dark.format(lat="collsc.bmad", root=root, sig_pz="8.804506566858e-08",
                           extra=load.format(f=beam_file) + sr_only + extra)
        (wd / f"{root}.nml").write_text(to_groups(text))
        r = subprocess.run([str(exe), f"{root}.nml"], cwd=wd, capture_output=True, text=True,
                           env={"OMP_NUM_THREADS": "4", "PATH": "/usr/bin:/bin"})
        if expect_fail:
            return r
        if r.returncode != 0:
            print(f"FAIL: {root} exited {r.returncode}:\n{r.stdout[-1500:]}")
            sys.exit(1)
        d = np.loadtxt(wd / f"{root}.diag.txt")
        ns = int(d[:, 1].max())
        return d.reshape(-1, ns, d.shape[1])

    collapse_slice(wd / "collsc_src-initial.beam.h5", wd / "collpc_one.beam.h5", 0, 1)
    collapse_slice(wd / "collsc_src-initial.beam.h5", wd / "collpc_two.beam.h5", 0, 2)
    collapse_slice(wd / "collsc_src-initial.beam.h5", wd / "collpc_nq.beam.h5", 0, 4,
                   zero_charge=True)

    # No radial scale to work from, and the solve says so rather than crashing or
    # returning a zero. The same beam with a scale stated runs.
    r = sc_run("collpc_r0", "collpc_one.beam.h5", "  space_charge%rmax = 0\n", expect_fail=True)
    refused = r.returncode != 0 and "NO RADIAL SCALE" in (r.stdout + r.stderr)
    r_ok = refused
    print(f"--- short-range with no radial scale: refused (exit {r.returncode})  "
          f"{'ok' if r_ok else 'FAIL'}")

    # The configuration domains, refused where they are stated.
    for tag, extra, phrase in (
            ("collpc_nphi", "  space_charge%rmax = 250e-6\n  space_charge%nphi = -1\n",
             "SPACE_CHARGE%NPHI"),
            ("collpc_ngrid", "  space_charge%rmax = 250e-6\n  space_charge%ngrid = 1\n",
             "SPACE_CHARGE%NGRID")):
        r = sc_run(tag, "collsc_src-initial.beam.h5", extra, expect_fail=True)
        got = r.returncode != 0 and phrase in (r.stdout + r.stderr)
        r_ok = r_ok and got
        print(f"--- {phrase.lower()} outside its domain: refused at setup  "
              f"{'ok' if got else 'FAIL'}")

    # Each case against itself with the solve off, so what is compared is this term's own
    # work and not the radiation the dark beam exchanges with its field either way.
    scale = "  space_charge%rmax = 250e-6\n"
    got = {}
    for tag, bf in (("one", "collpc_one.beam.h5"), ("two", "collpc_two.beam.h5"),
                    ("nq", "collpc_nq.beam.h5")):
        on = sc_run(f"collpc_{tag}", bf, scale)
        off = sc_run(f"collpc_{tag}f", bf,
                     scale + "  bmad_com%csr_and_space_charge_on = F\n")
        got[tag] = (float((on[-1, 0, 6] - on[0, 0, 6]) - (off[-1, 0, 6] - off[0, 0, 6])),
                    float((on[-1, 0, 7] - on[0, 0, 7]) - (off[-1, 0, 7] - off[0, 0, 7])))
    # The collapsed slice is slice 1. A lone particle feels its own field and nothing
    # else, and its harmonic sum cancels exactly; a slice with no charge has no source at
    # all. Neither moves in energy. Two colocated particles at different phases push each
    # other, equally and oppositely, so their mean holds while their spread opens.
    d_one, s_one = abs(got["one"][0]), abs(got["one"][1])
    d_nq, s_nq = abs(got["nq"][0]), abs(got["nq"][1])
    d_two, s_two = abs(got["two"][0]), abs(got["two"][1])
    # The pair's mean is not exactly still: the two exchange energy, which moves their
    # phases apart, so the push stops being exactly equal and opposite within the step.
    # Measured it is a percent of the spread they open, which is the statement that this
    # is a mutual force and not a common one.
    c_ok = (d_one < 1e-6 and s_one < 1e-6 and d_nq < 1e-6 and s_nq < 1e-6
            and s_two > 1e-3 and d_two < 5e-2 * s_two)
    ok = ok and r_ok and c_ok
    print(f"--- a collapsed slice on a grid with a scale, the solve's own work: lone "
          f"particle dE {d_one:.1e} eV and dsigma {s_one:.1e} eV, zero-charge slice dE "
          f"{d_nq:.1e} eV, colocated pair dsigma {s_two:.3e} eV against its mean "
          f"{d_two:.1e} eV  {'ok' if c_ok else 'FAIL'}")

    print("collective checks:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
