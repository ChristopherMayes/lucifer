#!/usr/bin/env python3
"""
Read a Lucifer particle dump into per-slice arrays, in either format.

The two formats carry the same beam in different charts, so the checks should not each
learn both. Genesis .par.h5 stores (x, y, px, py, theta, gamma) per slice with one
current, and openPMD .beam.h5 stores position, momentum, time and a per-particle
weight, with one particlePatch per slice and nothing else describing the window.

What comes back per slice: x, y [m], px, py [gamma*beta, Genesis's units, which
openPMD reaches through the file's own p0c], gamma, theta [rad], weight [C],
current [A], and n. Empty slices come back with n = 0 and empty arrays.

An openPMD dump adds pz, the Bmad energy deviation (p - p0)/p0, exactly as the file
states it. gamma is the same quantity through gamma = sqrt((p0(1+pz)/mc)^2 + 1), which
costs an ulp of gamma, and a check that differences energies between two runs pays that
ulp against a difference many orders smaller: a wake kick of 3e-8 in pz read back
through gamma agrees to 1.2e-8 rather than to rounding. Such a check should read pz. A
Genesis dump states gamma alone and carries no reference momentum, so pz is not
available there and the key is absent. Nor is it available from an openPMD file written
against the lattice's reference rather than its own, which is what the converter writes.

One caveat, about the wavelength. An openPMD dump does not carry it, since it belongs to
the run rather than to the beam. theta scales with it, so a caller reading an openPMD file
must say which wavelength its deck used. A Genesis file states its own and the argument is
then a cross-check.

theta is the absolute ponderomotive phase on both paths. Neither format has a place for a
run's reference phase phi0 and every reader restarts it at zero, so each writer folds it
into what it stores: Genesis stores theta itself, and an openPMD beam file stores the lag
theta implies, time = -theta/(ks c). The two therefore mean the same thing here, and a
dump is a restart point that reproduces the beam's phase against the field's.
"""

from __future__ import annotations

import pathlib

import h5py
import numpy as np

C_LIGHT = 2.99792458e8
M_ELECTRON = 0.51099895069e6      # eV, Bmad's value
E_CHARGE = 1.602176634e-19


def _record(group, name, n):
    """A dataset, or the constant a pseudo-dataset stands for, as an array of length n."""
    obj = group[name]
    if isinstance(obj, h5py.Dataset) and obj.shape:
        return np.asarray(obj[...], dtype=float)
    return np.full(n, float(np.atleast_1d(obj.attrs["value"])[0]))


def _attr(h5, name):
    return np.atleast_1d(h5.attrs[name])[0]


def is_openpmd(path):
    """True for an openPMD file, by the root attribute the standard requires."""
    with h5py.File(path) as h5:
        return "openPMD" in h5.attrs


def n_slice(path):
    """The window's slice count, in either format: slicecount, or the patch count."""
    with h5py.File(path) as h5:
        if "openPMD" not in h5.attrs:
            return int(h5["slicecount"][0])
        name = sorted(h5["data"].keys())[0]
        sp = h5[f"data/{name}/particles"]
        g = sp[sorted(sp.keys())[0]]
        return int(np.atleast_1d(g["particlePatches/numParticles"][...]).size)


def genesis_window(path):
    """
    The window a Genesis dump states: wavelength, slice spacing and slice count.

    An openPMD dump carries the count and nothing else, since the wavelength and the
    spacing belong to the run rather than to the beam. So when the two codes' dumps are
    compared, this is where the window comes from: the Genesis file states it, and the
    openPMD file beside it is read on those terms.
    """
    with h5py.File(path) as h5:
        if "openPMD" in h5.attrs:
            raise ValueError(f"{pathlib.Path(path).name} is openPMD and states no window: "
                             "read the window from the Genesis dump it is compared against")
        return dict(wavelength=float(h5["slicelength"][0]),
                    spacing=float(h5["slicespacing"][0]),
                    nslice=int(h5["slicecount"][0]))


def read_slices(path, wavelength=None, spacing=None):
    """Per-slice arrays from a Lucifer particle dump, either format.

    wavelength [m] is required for an openPMD dump and optional for a Genesis one, where
    a value that disagrees with the file is an error rather than a silent preference.
    """
    path = pathlib.Path(path)
    if is_guiding_centre(path):
        if wavelength is None:
            raise ValueError(f"{path.name} is a guiding-centre frame and does not carry the "
                             "wavelength: pass wavelength= (the deck's lambda0)")
        return _read_guiding_centre(path, wavelength, spacing)
    if is_openpmd(path):
        if wavelength is None:
            raise ValueError(f"{path.name} is openPMD and does not carry the wavelength: "
                             "pass wavelength= (the deck's lambda0)")
        return _read_openpmd(path, wavelength, spacing)
    return _read_genesis(path, wavelength)


def is_guiding_centre(path):
    """Whether the file is the diagnostic frame an averaged undulator's interior writes,
    the map's own chart in a file of its own (doc/reading-output.md)."""
    with h5py.File(path) as h5:
        return "guidingCentre" in h5


def _read_guiding_centre(path, wavelength, spacing=None):
    """The packed chart as the tracker held it, in the arrays the other readers return.
    The momenta come back in gamma*beta units through the file's own p0c, gamma from pz,
    and theta as phi0 + ks z / beta, which is the phase the tracker used with no fold in
    between. No conversion is applied: these are the coordinates the map evolved."""
    ks = 2 * np.pi / wavelength
    with h5py.File(path) as h5:
        g = h5["guidingCentre"]
        p0c = float(np.atleast_1d(g.attrs["p0c"])[0])
        phi0 = float(np.atleast_1d(g.attrs["phi0"])[0])
        counts = np.atleast_1d(g["sliceCount"][...]).astype(int)
        cols = {k: np.atleast_1d(g[k][...]) for k in ("x", "px", "y", "py", "z", "pz", "weight", "id")}
    p0_mc = p0c / M_ELECTRON
    if spacing is None:
        spacing = wavelength
    out, off = [], 0
    for n in counts:
        n = int(n)
        sl = slice(off, off + n)
        p_mc = p0_mc * (1 + cols["pz"][sl])
        gamma = np.sqrt(p_mc ** 2 + 1)
        beta = p_mc / gamma
        out.append(dict(n=n, x=cols["x"][sl], y=cols["y"][sl],
                        px=cols["px"][sl] * p0_mc, py=cols["py"][sl] * p0_mc,
                        gamma=gamma, theta=phi0 + ks * cols["z"][sl] / beta,
                        weight=cols["weight"][sl], pz=cols["pz"][sl], id=cols["id"][sl],
                        current=float(cols["weight"][sl].sum()) * C_LIGHT / spacing))
        off += n
    return out


def _read_genesis(path, wavelength=None):
    out = []
    with h5py.File(path) as h5:
        nslice = int(h5["slicecount"][0])
        spacing = float(h5["slicespacing"][0])
        lam = float(h5["slicelength"][0])
        if wavelength is not None and abs(wavelength - lam) > 1e-12 * lam:
            raise ValueError(f"{path.name} states wavelength {lam:.6e} m, caller said "
                             f"{wavelength:.6e} m")
        for i in range(1, nslice + 1):
            g = h5[f"slice{i:06d}"]
            n = g["gamma"].shape[0]
            current = float(g["current"][0])
            # The format has no per-particle weight: it carries one current, which the
            # reader divides out uniformly. That is the whole reason for .beam.h5.
            w = np.full(n, current * spacing / (C_LIGHT * n)) if n else np.zeros(0)
            out.append(dict(n=n, x=g["x"][...], y=g["y"][...],
                            px=g["px"][...], py=g["py"][...],
                            gamma=g["gamma"][...], theta=g["theta"][...],
                            weight=w, current=current))
    return out


def _read_openpmd(path, wavelength, spacing=None):
    """The particlePatches are the slices: one per slice, in window order, empty ones as
    zero-count patches. Nothing else in the file describes the window."""
    ks = 2 * np.pi / wavelength
    if spacing is None:
        spacing = wavelength
    out = []
    with h5py.File(path) as h5:
        name = sorted(h5["data"].keys())[0]
        sp = h5[f"data/{name}/particles"]
        g = sp[sorted(sp.keys())[0]]
        n_tot = int(np.atleast_1d(g.attrs["numParticles"])[0])
        pp = g["particlePatches"]
        n_pat = np.atleast_1d(pp["numParticles"][...]).astype(int)
        off = np.atleast_1d(pp["numParticlesOffset"][...]).astype(int)

        # Two shapes of the same beam. Bmad's writer states a reference momentum, so the
        # energy comes back as the deviation the tracker keeps. A file written against the
        # lattice's reference states none, which is deliberate (the converter leaves the
        # run's one reference to the lattice), and then only the total momentum is
        # available: gamma follows from it and pz has no reference to be measured against.
        pz = None
        if "totalMomentum" in g:
            p0c = _record(g, "totalMomentumOffset", n_tot)
            pz = _record(g, "totalMomentum", n_tot) / p0c    # written as pz * p0c
            p_mc = p0c * (1 + pz) / M_ELECTRON
        else:
            p_mc = np.sqrt(sum(_record(g, f"momentum/{c}", n_tot) ** 2
                               for c in ("x", "y", "z"))) / M_ELECTRON
        w = _record(g, "weight", n_tot)
        # time is t - t_ref = -theta/(ks c), the whole phase including the run's
        # reference. See the module docstring.
        theta = -ks * C_LIGHT * _record(g, "time", n_tot)
        x, y = _record(g, "position/x", n_tot), _record(g, "position/y", n_tot)
        px = _record(g, "momentum/x", n_tot) / M_ELECTRON
        py = _record(g, "momentum/y", n_tot) / M_ELECTRON
        gamma = np.sqrt(p_mc ** 2 + 1)

        for k, n in enumerate(n_pat):
            s = slice(off[k], off[k] + n)
            out.append(dict(n=int(n), x=x[s], y=y[s], px=px[s], py=py[s],
                            gamma=gamma[s], theta=theta[s], weight=w[s],
                            current=float(w[s].sum()) * C_LIGHT / spacing,
                            **({} if pz is None else dict(pz=pz[s]))))
    return out


FRAME_FORMAT = "lucifer-frames 2.0"


def unquiver(path, slices, wavelength):
    """A frame's slices with the undulator quiver taken back off, the guiding centre the
    averaged map tracks.

    The writer rebuilds the orbit into a .beam.h5 written inside an averaged undulator when
    a deck asks for the export, since an openPMD momentum record is the instantaneous
    kinetic momentum (fel-physics.md sec-quiverdump). This is that conversion's inverse, in
    the order the page states: the momenta first, which depend on the position in the
    device and the roll-off alone, then the longitudinal momentum from the recovered
    average, then the position and the lag.

    A frame that holds the guiding centre already comes back unchanged: a .gc.h5, and a
    frame the writer did not convert, one from a break, one this mode wrote with the orbit
    already in it, and one at an element face, where the device's field has ended.

    The frame's attributes are read from its group lucifer/frame, and an interior frame
    at any frameFormat but FRAME_FORMAT raises rather than being read under a layout this
    inverse does not know.
    """
    if is_guiding_centre(path):
        return slices
    with h5py.File(path) as h5:
        fr = h5.get("lucifer/frame")
        if fr is None or "aw" not in fr.attrs:
            return slices
        method = _attr(fr, "felMethod") if "felMethod" in fr.attrs else b""
        if isinstance(method, bytes):
            method = method.decode()
        if str(method) == "Unaveraged":
            return slices
        fmt = _attr(fr, "frameFormat") if "frameFormat" in fr.attrs else b""
        fmt = fmt.decode() if isinstance(fmt, bytes) else str(fmt)
        aw, ku = float(_attr(fr, "aw")), float(_attr(fr, "ku"))
        tilt, helical = float(_attr(fr, "tilt")), int(_attr(fr, "helical")) == 1
        s_ele, l_ele = float(_attr(fr, "sElement")), float(_attr(fr, "elementLength"))
        l_ramp = float(_attr(fr, "rampPeriods")) * 2 * np.pi / ku
    if s_ele <= 0 or s_ele >= l_ele:
        return slices
    if fmt != FRAME_FORMAT:
        raise ValueError(f"{pathlib.Path(path).name} is an interior frame of an averaged "
                         f"undulator at {fmt or 'no frameFormat'}, and this inverse reads "
                         f"{FRAME_FORMAT}")

    g, a_cos, a_sin, b_cos = _quiver_integrals(s_ele, l_ele, l_ramp, ku)
    a0 = aw if helical else np.sqrt(2.0) * aw
    ct, st = np.cos(tilt), np.sin(tilt)
    kx, ky = (ku * ku / 2, ku * ku / 2) if helical else (0.0, ku * ku)
    out = []
    for sl in slices:
        if sl["n"] == 0:
            out.append(sl)
            continue
        d = dict(sl)
        # The writer took the roll-off at the guiding centre, so the inverse iterates it:
        # the first pass reads it off the record's own position and the second off the
        # position that pass recovered. A planar device needs no second pass, its roll-off
        # depending on the wiggle frame's y, which the quiver does not move. A helical one
        # does, and the pass takes the inverse from second order in the quiver to fourth.
        xg, yg = sl["x"], sl["y"]
        for _ in range(2):
            xl, yl = ct * xg + st * yg, -st * xg + ct * yg
            awloc = 1 + 0.5 * (kx * xl * xl + ky * yl * yl)
            amp = a0 * awloc
            ux, uy = ct * sl["px"] + st * sl["py"], -st * sl["px"] + ct * sl["py"]
            ux = ux - amp * g * np.cos(ku * s_ele)
            if helical:
                uy = uy - amp * g * np.sin(ku * s_ele)
            us = np.sqrt(sl["gamma"] ** 2 - 1 - ux * ux - uy * uy - (aw * awloc) ** 2)
            dxl = amp / us * a_cos
            dyl = amp / us * a_sin if helical else np.zeros_like(dxl)
            xg = sl["x"] - (ct * dxl - st * dyl)
            yg = sl["y"] - (st * dxl + ct * dyl)
        if helical:
            dtau = sl["gamma"] / (2 * us ** 3) * 2 * amp * (ux * a_cos + uy * a_sin)
        else:
            dtau = sl["gamma"] / (2 * us ** 3) * (2 * ux * amp * a_cos + (aw * awloc) ** 2 * b_cos)
        d["x"], d["y"] = xg, yg
        d["px"] = ct * ux - st * uy
        d["py"] = st * ux + ct * uy
        # theta = phi0 + ks z / beta and the writer took beta * dtau off z, so the lag
        # comes off the phase as ks * dtau, with no beta left in it.
        d["theta"] = sl["theta"] + 2 * np.pi / wavelength * dtau
        out.append(d)
    return out


def _quiver_envelope(t, l_ele, l_ramp):
    """The device's amplitude envelope, the reader's side of fel_und_envelope."""
    g = np.ones_like(t)
    if l_ramp > 0:
        m = t < l_ramp;            g[m] = np.sin(np.pi * t[m] / (2 * l_ramp)) ** 2
        m = t > l_ele - l_ramp;    g[m] = np.sin(np.pi * (l_ele - t[m]) / (2 * l_ramp)) ** 2
    return g


def _quiver_integrals(s_ele, l_ele, l_ramp, ku):
    """g at s and the three phase integrals from the upstream face, Simpson at the step the
    writer uses. The ramp is what makes them worth integrating rather than expanding."""
    lam_u = 2 * np.pi / ku
    n = max(64, 2 * int(round(128 * s_ele / lam_u)))
    t = np.linspace(0.0, s_ele, n + 1)
    w = np.ones(n + 1); w[1:-1:2] = 4; w[2:-1:2] = 2
    g = _quiver_envelope(t, l_ele, l_ramp)
    h = s_ele / n
    return (float(_quiver_envelope(np.array([s_ele]), l_ele, l_ramp)[0]),
            h / 3 * float(np.dot(w, g * np.cos(ku * t))),
            h / 3 * float(np.dot(w, g * np.sin(ku * t))),
            h / 3 * float(np.dot(w, g * g * np.cos(2 * ku * t))))
