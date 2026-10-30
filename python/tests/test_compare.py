"""lucifer.compare without the executable: a Genesis4 output of the file's layout, built
here, read into Curves, differenced against a run, and refused where the two runs do not
line up. The figure is drawn where matplotlib is installed.

    python -m pytest python/tests -q
"""
import pathlib

import h5py
import numpy as np
import pytest

from lucifer import compare as cmp

M_ELECTRON = 0.51099895069e6
NZ, NS = 12, 4


def curves(name="run", nz=NZ, ns=NS, scale=1.0, ds=3e-10, kind="lucifer"):
    """Curves of a made-up run: power growing along s, one slice dark at the start."""
    s = np.linspace(0.0, 4.0, nz)
    power = np.outer(np.exp(s), 1.0 + np.arange(ns)) * 1e3 * scale
    power[0, 0] = 0.0
    bunching = np.outer(s / 4.0, np.ones(ns)) * 0.1 * scale
    energy = 5.8e9 - np.outer(s ** 2, np.ones(ns)) * 1e3 * scale
    size = np.full((nz, ns), 2e-5)
    return cmp.Curves(name, kind, s, power, bunching, energy, size, size * 1.5,
                      ds if ns > 1 else None)


def genesis_out(path, c):
    """A Genesis4 .out.h5 holding the curves ``c``, energy as gamma."""
    with h5py.File(path, "w") as h5:
        h5["Lattice/zplot"] = c.s
        h5["Field/power"] = c.power
        h5["Beam/bunching"] = c.bunching
        h5["Beam/energy"] = c.energy / M_ELECTRON
        h5["Beam/xsize"] = c.size_x
        h5["Beam/ysize"] = c.size_y
        h5["Global/s"] = np.arange(c.n_slice) * (c.ds or 0.0)


def genesis_field(path, u, dx):
    """A Genesis4 .fld.h5 of the field u (nslice, n, n), in the dfl scale."""
    from lucifer.field import MU0_C
    dfl = u * dx / np.sqrt(2 * MU0_C)
    with h5py.File(path, "w") as h5:
        h5["gridpoints"] = np.array([u.shape[1]], dtype=np.int32)
        h5["gridsize"] = np.array([dx])
        h5["refposition"] = np.array([0.0])
        h5["wavelength"] = np.array([1e-10])
        h5["slicecount"] = np.array([u.shape[0]], dtype=np.int32)
        h5["slicespacing"] = np.array([3e-10])
        for i in range(u.shape[0]):
            g = h5.create_group(f"slice{i + 1:06d}")
            g["field-real"] = np.ascontiguousarray(dfl[i].real).reshape(-1)
            g["field-imag"] = np.ascontiguousarray(dfl[i].imag).reshape(-1)


class FakeStats:
    """What from_stats reads of a Stats, for a run with no file."""

    def __init__(self, c, path):
        self.c, self.path = c, pathlib.Path(path)
        self.s = c.s
        self.ct_slice = np.arange(c.n_slice) * c.ds

    def __getitem__(self, key):
        nz, ns = self.c.power.shape
        if key == "beam/slice/sigma":
            sig = np.zeros((nz, ns, 6, 6))
            sig[:, :, 0, 0] = self.c.size_x ** 2
            sig[:, :, 2, 2] = self.c.size_y ** 2
            return sig
        return {"field/total/power": self.c.power, "beam/slice/bunching": self.c.bunching,
                "beam/slice/energy": self.c.energy}[key]


def test_a_genesis_output_is_read_in_ev_with_its_slice_spacing(tmp_path):
    c = curves()
    genesis_out(tmp_path / "Aramis.out.h5", c)
    g = cmp.read_curves(tmp_path / "Aramis.out.h5")
    assert g.kind == "genesis4" and g.name == "Genesis4"
    assert np.allclose(g.energy, c.energy) and g.ds == pytest.approx(3e-10)
    assert g.pulse_energy is not None and g.pulse_energy.shape == (NZ,)
    assert repr(g) == "Curves(Genesis4, 12 records to s = 4 m, 4 slices)"


def test_a_stats_file_gives_its_sizes_from_sigma_and_its_name_from_its_stem():
    c = curves()
    r = cmp.from_stats(FakeStats(c, "s7.stats.h5"))
    assert r.name == "s7" and r.kind == "lucifer"
    assert np.allclose(r.size_x, c.size_x) and np.allclose(r.size_y, c.size_y)
    assert r.ds == pytest.approx(3e-10)


def test_a_file_of_neither_kind_is_refused_with_its_roots(tmp_path):
    with h5py.File(tmp_path / "other.h5", "w") as h5:
        h5["data/x"] = np.zeros(3)
    with pytest.raises(ValueError, match="neither a bmad-stats file nor a Genesis4 output: its root holds data"):
        cmp.read_curves(tmp_path / "other.h5")


def test_aligned_names_what_differs():
    a = curves()
    assert cmp.aligned(a, curves("ref")) is None
    assert "records" in cmp.aligned(a, curves("ref", nz=NZ + 1))
    b = curves("ref")
    b.s = b.s + 1e-6
    assert "differ by up to 1e-06 m" in cmp.aligned(a, b)
    assert "slices" in cmp.aligned(a, curves("ref", ns=NS - 1))
    assert "slice spacings" in cmp.aligned(a, curves("ref", ds=4e-10))
    with pytest.raises(ValueError, match="slices"):
        cmp.differences(a, curves("ref", ns=NS - 1))


def test_differences_are_divided_as_the_labels_say_and_reduced_over_the_slices():
    ref = curves("ref")
    run = curves("run", scale=1.0 + 1e-5)
    d = cmp.differences(run, ref)
    assert set(d) == {"power", "bunching", "energy", "exit_power"}
    # The power is elementwise relative, so a uniform factor is that factor everywhere
    # light is, and the dark first slice divides by the floor rather than by zero.
    assert d["power"].largest == pytest.approx(1e-5, rel=1e-6)
    assert np.isfinite(d["power"].values).all()
    # The bunching is over the reference's peak, so the factor is seen at the peak alone.
    assert d["bunching"].largest == pytest.approx(1e-5, rel=1e-6)
    assert d["bunching"].values[0] == 0.0
    # The energy change is over the reference's largest change.
    assert d["energy"].largest == pytest.approx(1e-5, rel=1e-6)
    # The exit power across the window, at the slice positions.
    assert d["exit_power"].x.shape == (NS,) and d["exit_power"].x[1] == pytest.approx(3e-10)
    assert d["exit_power"].largest == pytest.approx(1e-5, rel=1e-6)
    assert repr(d["power"]).startswith("Difference(power, largest 1.000e-05")


def test_a_steady_state_run_differences_its_sizes():
    ref, run = curves("ref", ns=1), curves("run", ns=1)
    run.size_x = run.size_x * 1.01
    d = cmp.differences(run, ref)
    assert set(d) == {"power", "bunching", "energy", "size_x", "size_y"}
    assert d["size_x"].largest == pytest.approx(0.01) and d["size_y"].largest == 0.0
    assert d["power"].values.shape == (NZ,)


def test_the_sign_of_the_largest_slice_difference_is_kept():
    ref = curves("ref")
    run = curves("run")
    run.bunching[5, 2] = ref.bunching[5, 2] - 0.3 * np.abs(ref.bunching).max()
    run.bunching[5, 1] = ref.bunching[5, 1] + 0.1 * np.abs(ref.bunching).max()
    d = cmp.differences(run, ref)["bunching"]
    assert d.values[5] == pytest.approx(-0.3)


def test_the_field_row_and_level_come_from_the_dumps(tmp_path):
    n = 9
    x = np.linspace(-2, 2, n)
    u = (np.exp(-x[:, None] ** 2 - x[None, :] ** 2) * (1 + 0.5j))[None]
    genesis_field(tmp_path / "ref.fld.h5", u, 1e-6)
    genesis_field(tmp_path / "run.fld.h5", u * (1 + 2e-4), 1e-6)
    ref, run = curves("ref", ns=1), curves("run", ns=1)
    fields = (tmp_path / "run.fld.h5", tmp_path / "ref.fld.h5")
    d = cmp.differences(run, ref, fields=fields)["field"]
    assert d.x.shape == (n,) and d.x[n // 2] == 0.0 and d.x[1] - d.x[0] == pytest.approx(1e-6)
    assert d.largest == pytest.approx(2e-4, rel=1e-6)
    assert cmp.field_level(fields) == pytest.approx(2e-4, rel=1e-6)
    genesis_field(tmp_path / "coarse.fld.h5", u[:, ::2, ::2], 2e-6)
    with pytest.raises(ValueError, match="field grids differ"):
        cmp.differences(run, ref, fields=(tmp_path / "coarse.fld.h5", tmp_path / "ref.fld.h5"))


def test_compare_draws_four_lettered_blocks_with_a_power_of_ten_on_each_difference(capsys):
    pytest.importorskip("matplotlib")
    pytest.importorskip("beamphysics")
    from lucifer.plot import compare
    ref = curves("Genesis4", kind="genesis4")
    run = curves("s7", scale=1.0 + 1e-5)
    fig = compare(run, ref, return_figure=True)
    letters = sorted(t.get_text() for a in fig.axes for t in a.texts if t.get_text().startswith("("))
    assert letters == ["(a)", "(b)", "(c)", "(d)"]
    labels = [a.get_ylabel() for a in fig.axes]
    assert any(r"\delta P / P" in lab and "10^{-5}" in lab for lab in labels), labels
    assert any(lab.startswith("$E_\\mathrm{pulse}$") for lab in labels), labels
    # A run over Genesis4 is Lucifer in the legend, whatever its file was called.
    assert [t.get_text() for t in fig.axes[0].get_legend().get_texts()] == ["Genesis4", "Lucifer"]
    # Misaligned runs draw the overlays alone and say why.
    short = curves("ref", nz=NZ + 1)
    fig2 = compare(run, short, return_figure=True)
    assert len(fig2.axes) == len(fig.axes) - 4
    assert "the differences are not drawn" in capsys.readouterr().out
