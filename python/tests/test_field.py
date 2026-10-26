"""lucifer.field without the executable: an escaped file of the program's layout, built
here, read back whole, and a file without the private group refused.

    python -m pytest python/tests -q
"""
import h5py
import numpy as np
import pytest

from lucifer.field import ESCAPED_GROUP, PARAM_NAMES, read_escaped

NX, NZ = 8, 5


def escaped_file(path, two_pol=False):
    """An escaped file as the program writes it: an EXT_Wavefront mesh with a NaN
    zCoordinate, and lucifer/escaped beside it."""
    from beamphysics.wavefront import Wavefront
    x = np.linspace(-2, 2, NX)
    env = np.exp(-x[:, None, None] ** 2 - x[None, :, None] ** 2) * (1 + np.arange(NZ))[None, None, :]
    ey = env * 0.5 if two_pol else None
    Wavefront(Ex=env.astype(complex), Ey=None if ey is None else ey.astype(complex),
              dx=1e-5, dy=1e-5, dz=1e-10, wavelength=1e-10, s_position=float("nan")).write_openpmd(str(path))
    z = np.linspace(1.3, 0.1, NZ)
    pms = np.arange(NZ * 25, dtype=float).reshape(NZ, 25)
    with h5py.File(path, "a") as h5:
        g = h5.create_group(ESCAPED_GROUP)
        g["zTransmit"] = z
        g["wavefrontParams"] = pms
        g["wavefrontParams"].attrs["columns"] = np.array(PARAM_NAMES, dtype="S")
        if two_pol:
            g["wavefrontParamsY"] = pms * 2
            g["wavefrontParamsY"].attrs["columns"] = np.array(PARAM_NAMES, dtype="S")
    return z, pms


def test_read_escaped_returns_the_planes_and_the_two_arrays(tmp_path):
    pytest.importorskip("beamphysics")
    z, pms = escaped_file(tmp_path / "s7-escaped.wf.h5")
    esc = read_escaped(tmp_path / "s7-escaped.wf.h5")
    assert esc.wavefront.nz == NZ and esc.wavefront.Ex.shape == (NX, NX, NZ)
    assert np.isnan(esc.wavefront.s_position)
    assert np.array_equal(esc.z_transmit, z) and np.array_equal(esc.params, pms)
    assert esc.params_y is None
    assert repr(esc).startswith(f"Escaped({NZ} planes, left at 0.1 to 1.3 m")


def test_a_second_polarization_brings_its_params(tmp_path):
    pytest.importorskip("beamphysics")
    _, pms = escaped_file(tmp_path / "s7-escaped.wf.h5", two_pol=True)
    esc = read_escaped(tmp_path / "s7-escaped.wf.h5")
    assert esc.wavefront.Ey is not None and np.array_equal(esc.params_y, pms * 2)


def test_a_file_without_the_group_is_refused(tmp_path):
    pytest.importorskip("beamphysics")
    escaped_file(tmp_path / "s7-escaped.wf.h5")
    with h5py.File(tmp_path / "s7-escaped.wf.h5", "a") as h5:
        del h5["lucifer"]
    with pytest.raises(ValueError, match=ESCAPED_GROUP):
        read_escaped(tmp_path / "s7-escaped.wf.h5")


def test_other_columns_are_refused(tmp_path):
    pytest.importorskip("beamphysics")
    escaped_file(tmp_path / "s7-escaped.wf.h5")
    with h5py.File(tmp_path / "s7-escaped.wf.h5", "a") as h5:
        h5[ESCAPED_GROUP + "/wavefrontParams"].attrs["columns"] = np.array(["other"] * 25, dtype="S")
    with pytest.raises(ValueError, match="names its params columns"):
        read_escaped(tmp_path / "s7-escaped.wf.h5")
