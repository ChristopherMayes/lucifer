"""The dumps of a run without the executable: the keys the stats file's list names, a frame series
expanded, the dumps read on first access and kept, one slice read from its patch, and an
archive carrying its dumps. Every file is written here.

    python -m pytest python/tests -q
"""
import pathlib

import h5py
import numpy as np
import pytest

from lucifer import Lucifer
from lucifer.dumps import Dumps, dump_files, expand
from lucifer.runner import Result

FILES = ["s7-initial.beam.h5", "s7-initial.wf.h5", "s7-at1-UND.beam.h5", "s7-at1-UND.wf.h5",
         "s7-final.beam.h5", "s7-final.wf.h5", "s7-final-h3.wf.h5", "s7-escaped.wf.h5", "s7-pulse.wf.h5",
         "s7-{000001..000003}.beam.h5", "s7-{000001..000003}.wf.h5", "s7-{000001..000003}-h3.wf.h5",
         "s7.stats.h5", "s7.diag.txt"]
COUNTS = [3, 5, 2]


def beam_file(path):
    """A beam dump of three slices, the program's layout: one species with its particle
    patches, numParticles and numParticlesOffset."""
    from beamphysics import ParticleGroup
    n = sum(COUNTS)
    rng = np.random.default_rng(7)
    data = {k: rng.normal(size=n) for k in ("x", "px", "y", "py", "z")}
    data.update(pz=1e9 + rng.normal(size=n), t=rng.normal(size=n) * 1e-15,
                status=np.ones(n, dtype=int), weight=np.full(n, 1e-18), species="electron", id=np.arange(1, n + 1))
    ParticleGroup(data=data).write(str(path))
    from beamphysics.readers import particle_paths
    with h5py.File(path, "a") as h5:
        g = h5[particle_paths(h5)[0]]["electron"]
        g["particlePatches/numParticles"] = np.array(COUNTS)
        g["particlePatches/numParticlesOffset"] = np.concatenate([[0], np.cumsum(COUNTS)[:-1]])
    return path


def field_file(path):
    from beamphysics.wavefront import Wavefront
    Wavefront(Ex=np.ones((4, 4, 2), dtype=complex), dx=1e-5, dy=1e-5, dz=1e-10, wavelength=1e-10,
              s_position=1.0).write_openpmd(str(path))
    return path


def test_a_frame_series_expands_to_its_files():
    assert expand("s7-{000001..000003}.beam.h5") == ["s7-000001.beam.h5", "s7-000002.beam.h5", "s7-000003.beam.h5"]
    assert expand("s7-final.beam.h5") == ["s7-final.beam.h5"]


def test_the_stats_file_names_the_keys(tmp_path):
    # The list is the attribute files of run/, as the program writes it, and every file
    # it names is there.
    from lucifer.dumps import expand
    with h5py.File(tmp_path / "s7.stats.h5", "w") as f:
        f.create_group("run").attrs["files"] = np.array(FILES, dtype="S")
    for row in FILES:
        for name in expand(row):
            if name != "s7.stats.h5":
                (tmp_path / name).write_bytes(b"0")
    R = Result(0, "", tmp_path / "s7.toml", tmp_path, "s7")
    assert list(R.particle_files) == ["initial", "at1-UND", "final", 1, 2, 3]
    assert list(R.field_files) == ["initial", "at1-UND", "final", "final-h3", "escaped", "pulse",
                                   1, 2, 3, "1-h3", "2-h3", "3-h3"]
    assert R.particle_files[2] == tmp_path / "s7-000002.beam.h5"
    assert len(R.files) == 2 + 2 + 3 + 2 + 3 * 3 + 2


def test_an_out_root_with_dashes_keys_by_what_follows_it(tmp_path):
    beams, fields = dump_files([tmp_path / "a-b-final.beam.h5", tmp_path / "a-b-000012.wf.h5",
                                tmp_path / "a-b.stats.h5", tmp_path / "other-final.beam.h5"], "a-b")
    assert list(beams) == ["final"] and list(fields) == [12]


def test_a_dump_is_read_on_first_access_and_kept(tmp_path):
    reads = []
    # The loader is a stand-in, so the file need only exist: the test runs without
    # openpmd-beamphysics, as the docs workflow runs it.
    (tmp_path / "s7-final.wf.h5").write_bytes(b"0" * 2048)
    files = {"final": tmp_path / "s7-final.wf.h5"}
    d = Dumps(files, lambda p: reads.append(p) or p.name, "field")
    text = repr(d)
    assert "s7-final.wf.h5" in text and "on disk" in text and reads == []
    assert d["final"] == "s7-final.wf.h5" and d["final"] == "s7-final.wf.h5" and len(reads) == 1
    assert d.loaded == ["final"] and "loaded" in repr(d)
    with pytest.raises(KeyError, match="no field dump 'initial'"):
        d["initial"]


def test_the_loaders_read_the_files(tmp_path):
    pytest.importorskip("beamphysics")
    beam_file(tmp_path / "s7-final.beam.h5")
    field_file(tmp_path / "s7-final.wf.h5")
    R = Result(0, "", tmp_path / "s7.toml", tmp_path, "s7")
    assert len(R.particles["final"]) == sum(COUNTS)
    assert R.fields["final"].Ex.shape == (4, 4, 2)
    assert R.load_particles() == ["final"] and R.load_fields() == ["final"]


def test_one_slice_is_read_from_its_patch(tmp_path):
    pytest.importorskip("beamphysics")
    beam_file(tmp_path / "s7-final.beam.h5")
    R = Result(0, "", tmp_path / "s7.toml", tmp_path, "s7")
    whole = R.particles["final"]
    one = R.load_particles_by_key("final", slice=1)
    assert len(one) == COUNTS[1]
    for k in ("x", "px", "y", "py", "z", "pz", "t", "weight", "status"):
        assert np.array_equal(getattr(one, k), getattr(whole, k)[COUNTS[0]:COUNTS[0] + COUNTS[1]])
    with pytest.raises(IndexError, match="holds 3 slices"):
        R.load_particles_by_key("final", slice=3)


def test_an_archive_carries_its_dumps(tmp_path):
    pytest.importorskip("beamphysics")
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    with h5py.File(run_dir / "s7.stats.h5", "w") as f:
        f["meta/input_echo"] = b'lattice = "aramis.bmad"\n\n[run]\nout_root = "s7"\n'
    beam_file(run_dir / "s7-final.beam.h5")
    field_file(run_dir / "s7-final.wf.h5")
    L = Lucifer.from_archive(run_dir / "s7.stats.h5")
    L.output.dump_dir = run_dir
    bare = L.archive(tmp_path / "bare.h5")
    full = L.archive(tmp_path / "full.h5", include_dumps=True)
    assert sorted(p.name for p in (tmp_path / "full").iterdir()) == ["s7-final.beam.h5", "s7-final.wf.h5"]
    assert not (tmp_path / "bare").exists()
    assert list(Lucifer.from_archive(bare).output.particle_files) == []
    R = Lucifer.from_archive(full).output
    assert list(R.particle_files) == ["final"] and list(R.field_files) == ["final"]
    assert len(R.particles["final"]) == sum(COUNTS)


def test_a_temporary_directory_lives_as_long_as_the_result(tmp_path):
    import gc
    import tempfile
    d = tempfile.TemporaryDirectory(dir=tmp_path)
    where = pathlib.Path(d.name)
    R = Result(0, "", where / "s7.toml", where, "s7", keep=d)
    del d
    gc.collect()
    assert where.is_dir()
    del R
    gc.collect()
    assert not where.exists()
