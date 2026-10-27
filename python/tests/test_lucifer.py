"""The Lucifer object without the executable: construction, the input file it writes, its
identity, copying, restoring from an archive, its representations, and its names against
lume-base's CommandWrapper where lume-base is installed.

    python -m pytest python/tests -q
"""
import inspect
import pathlib
import tomllib

import pytest

from lucifer import Input, Lucifer

TEXT = 'lattice = "aramis.bmad"\n\n[beam]\nn_particle = 4096\n\n[run]\nout_root = "s7"\n'

# lume-base's names Lucifer leaves out: MPI, since the program threads with OpenMP, and
# from_yaml, since the input is TOML.
LEFT_OUT = {"use_mpi", "mpi_run", "command_mpi", "COMMAND_MPI", "MPI_RUN", "from_yaml"}


def test_version_is_the_installed_packages():
    # setuptools-scm sets the installed package's version from the git tag, so an install
    # carries a version and the package reports the one it was installed with.
    import importlib.metadata

    import lucifer
    assert isinstance(lucifer.__version__, str) and lucifer.__version__
    assert lucifer.__version__ == importlib.metadata.version("lucifer")


def test_built_from_an_input_a_path_and_text(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "s7.toml").write_text(TEXT)
    from_input = Lucifer(Input.loads(TEXT))
    from_path = Lucifer("s7.toml")
    from_text = Lucifer(TEXT)
    assert from_input.input == from_path.input == from_text.input
    assert from_path.input_file == str(tmp_path / "s7.toml")
    assert from_input.input_file is None


def test_the_written_input_names_its_files_by_absolute_paths(tmp_path):
    src = tmp_path / "case"
    src.mkdir()
    (src / "s7.toml").write_text(TEXT.replace('[beam]\n', '[beam]\nposition_file = "bunch.h5"\n'))
    L = Lucifer(src / "s7.toml", workdir=tmp_path / "runs")
    L.configure()
    written = pathlib.Path(L.input_file)
    assert written.parent.parent == tmp_path / "runs" and written.name == "s7.toml"
    data = tomllib.loads(written.read_text())
    assert data["lattice"] == str(src / "aramis.bmad")
    assert data["beam"]["position_file"] == str(src / "bunch.h5")
    assert L.configured and not L.finished


def test_the_input_files_own_directory_without_a_temporary_one(tmp_path):
    (tmp_path / "s7.toml").write_text(TEXT)
    L = Lucifer(tmp_path / "s7.toml", use_temp_dir=False)
    L.configure()
    assert L.path == str(tmp_path)


def test_fingerprint_copy_and_reset(tmp_path):
    L = Lucifer(TEXT, workdir=tmp_path)
    same = Lucifer(TEXT, workdir=tmp_path)
    assert L.fingerprint() == same.fingerprint()
    same.input.beam.n_particle = 2048
    assert L.fingerprint() != same.fingerprint()
    L.configure()
    C = L.copy()
    C.input.run.out_root = "other"
    assert L.input.run.out_root == "s7" and C.path is None and not C.configured
    L.reset()
    assert L.path is None and not L.configured


def test_restored_from_an_archive(tmp_path):
    h5py = pytest.importorskip("h5py")
    path = tmp_path / "run.h5"
    with h5py.File(path, "w") as f:
        f["meta/input_echo"] = Input.loads(TEXT).toml(explicit=True).encode()
    R = Lucifer.from_archive(path)
    assert R.input.beam.n_particle == 4096 and R.finished and R.output.stats_path == path
    assert repr(R) == "Lucifer(archive='run.h5', finished)"


def test_the_representations():
    inp = Input.loads(TEXT)
    assert repr(inp) == "Input(lattice='aramis.bmad', run.out_root='s7', beam.n_particle=4096)"
    assert "<table>" in inp._repr_html_() and "beam.n_particle" in inp._repr_html_()
    L = Lucifer(TEXT)
    assert repr(L) == "Lucifer(input='in memory', path=None, not run)"
    assert "not configured" in L._repr_html_()


def test_names_against_lume_base():
    """Every public name of lume-base's CommandWrapper is on Lucifer, less the MPI names and
    from_yaml. Skipped where lume-base is not installed, which is the usual case."""
    base = pytest.importorskip("lume.base", reason="lume-base is not installed, so the names are not compared")
    theirs = {n for n in dir(base.CommandWrapper) if not n.startswith("_")}
    ours = {n for n in dir(Lucifer) if not n.startswith("_")}
    missing = sorted(theirs - ours - LEFT_OUT)
    assert not missing, f"Lucifer lacks CommandWrapper's {missing}"
    for name in theirs & ours:
        a, b = inspect.getattr_static(base.CommandWrapper, name), inspect.getattr_static(Lucifer, name)
        assert isinstance(a, property) == isinstance(b, property), f"{name} is a property in one and not the other"
