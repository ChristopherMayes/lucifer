"""One object to build, run, read and plot a Lucifer run.

    from lucifer import Lucifer

    L = Lucifer("lucifer.toml", threads=12)
    L.input.beam.n_particle = 4096
    L.run()
    L.plot("field/total/power", yscale="log")
    L.archive("run.h5")
    L2 = Lucifer.from_archive("run.h5")

Lucifer carries the names and the meanings of lume-base's CommandWrapper, the interface
LUME-Genesis's Genesis4 and the other LUME codes share, without depending on lume-base:
input, output, initial_particles, workdir, path, command, timeout, verbose, configured,
finished, error, configure(), run(), write_input(), load_input(), load_output(),
fingerprint(), copy(), reset(), archive(), load_archive(), from_archive(), plot(). The
MPI names, use_mpi, mpi_run, command_mpi, COMMAND_MPI and MPI_RUN, are left out, since
the program threads with OpenMP, and from_yaml is left out, since the input is TOML.

A run happens in a directory of its own: a temporary directory inside ``workdir``, or
``workdir`` itself with ``use_temp_dir=False``, or the input file's own directory when
neither is given. The input file is written there with the lattice and every file it
names as absolute paths, resolved from the directory the input came from, so a run in a
temporary directory reads its files where they are and nothing is copied.

The archive of a run is its stats file, which holds the input file the run read in
meta/input_echo and every output array, so archive() copies it and from_archive()
rebuilds the object from it.
"""
from __future__ import annotations

import copy as _copy
import dataclasses
import hashlib
import os
import pathlib
import shutil
import tempfile

__all__ = ["Lucifer"]

# The keys that name a file, resolved to absolute paths when the input file is written
# into the run's directory.
PATH_KEYS = ("lattice", "beam.file", "beam.position_file")


class Lucifer:
    """A Lucifer run: its input, the directory it runs in, and its output.

    ``input`` is an Input, the path of an input file, or TOML text. ``workdir``,
    ``use_temp_dir``, ``command``, ``verbose`` and ``timeout`` mean what lume-base's
    CommandWrapper means by them, and ``threads`` sets OMP_NUM_THREADS. ``command`` is the
    executable, by default ``$LUCIFER_EXE`` or ``lucifer`` on PATH. ``initial_particles``
    is an openPMD-beamphysics ParticleGroup to start from."""

    COMMAND = "lucifer"
    WORKDIR = None

    def __init__(self, input=None, *, workdir=None, use_temp_dir=True, command=None, threads=None,
                 verbose=False, timeout=None, initial_particles=None):
        self._input = None
        self.original_path = None
        self.original_input_file = None
        self._output = None
        self._initial_particles = initial_particles
        self._command = command
        self.threads = threads
        self._verbose = verbose
        self._timeout = timeout
        self._use_temp_dir = use_temp_dir
        self._workdir = str(pathlib.Path(workdir).resolve()) if workdir else self.WORKDIR
        self._tempdir = None
        self._base_path = None
        self._archive = None
        self._configured = False
        self._finished = False
        self._error = False
        if input is not None:
            self._set_input(input)

    def _set_input(self, value):
        from lucifer.models import Input
        if isinstance(value, Input):
            self._input = value
            self.original_path = os.getcwd()
        elif isinstance(value, pathlib.Path) or ("\n" not in str(value) and pathlib.Path(value).is_file()):
            self.load_input(value)
        elif isinstance(value, str):
            self._input = Input.loads(value)
            self.original_path = os.getcwd()
        else:
            raise TypeError(f"input is an Input, a path or TOML text, not {type(value).__name__}")

    # ------------------------------------------------------------------
    # The state, CommandWrapper's properties.

    @property
    def input(self):
        """The input file, an Input."""
        return self._input

    @input.setter
    def input(self, value):
        self._set_input(value)

    @property
    def output(self):
        """The run's Result, with its stats file, or None before a run."""
        return self._output

    @output.setter
    def output(self, value):
        self._output = value

    @property
    def initial_particles(self):
        """An openPMD-beamphysics ParticleGroup the run starts from, or None."""
        return self._initial_particles

    @initial_particles.setter
    def initial_particles(self, value):
        self._initial_particles = value
        self._configured = False

    @property
    def command(self):
        """The executable: as given, else ``$LUCIFER_EXE``, else ``lucifer`` on PATH."""
        return self._command

    @command.setter
    def command(self, value):
        self._command = value

    @property
    def verbose(self):
        """Whether run() prints the program's output when it ends."""
        return self._verbose

    @verbose.setter
    def verbose(self, value):
        self._verbose = value

    @property
    def timeout(self):
        """Seconds after which a run is stopped and refused, or None."""
        return self._timeout

    @timeout.setter
    def timeout(self, value):
        self._timeout = value

    @property
    def configured(self):
        """Whether the run's directory exists and holds the input file."""
        return self._configured

    @configured.setter
    def configured(self, value):
        self._configured = value

    @property
    def finished(self):
        """Whether a run has ended, with or without an error."""
        return self._finished

    @finished.setter
    def finished(self, value):
        self._finished = value

    @property
    def error(self):
        """Whether the last run ended with an error."""
        return self._error

    @error.setter
    def error(self, value):
        self._error = value

    @property
    def use_temp_dir(self):
        """Whether the run happens in a temporary directory inside ``workdir``."""
        return self._use_temp_dir

    @property
    def workdir(self):
        """The directory runs happen in or under. Setting it sets up the run's directory."""
        return self._workdir

    @workdir.setter
    def workdir(self, value):
        self._workdir = str(pathlib.Path(value).resolve())
        self.setup_workdir(self._workdir)

    @property
    def path(self):
        """The run's directory, once configured."""
        return self._base_path

    @path.setter
    def path(self, value):
        self._base_path = value

    @property
    def input_file(self):
        """The input file the run reads: the one written into the run's directory once
        configured, else the file the input was loaded from."""
        if self._configured and self._base_path:
            return str(pathlib.Path(self._base_path) / self._input_name())
        if self.original_path and self.original_input_file:
            return str(pathlib.Path(self.original_path) / self.original_input_file)
        return None

    def _input_name(self):
        return self.original_input_file or "lucifer.toml"

    # ------------------------------------------------------------------
    # Setting up and running.

    @staticmethod
    def input_parser(path):
        """The Input of an input file."""
        from lucifer.models import Input
        return Input.read(path)

    def load_input(self, input_filepath, **kwargs):
        """Read an input file into ``input``, and remember its directory, from which the
        files it names are resolved."""
        f = pathlib.Path(input_filepath).resolve()
        self.original_path, self.original_input_file = str(f.parent), f.name
        self._input = self.input_parser(f)
        self._configured = False

    def setup_workdir(self, workdir, cleanup=True):
        """Set the run's directory: a temporary directory inside ``workdir`` with
        ``use_temp_dir``, else ``workdir``, else the input file's own directory."""
        if self._use_temp_dir:
            if workdir:
                pathlib.Path(workdir).mkdir(parents=True, exist_ok=True)
            self._tempdir = tempfile.TemporaryDirectory(dir=workdir, prefix="lucifer_")
            self._base_path = self._tempdir.name
        elif workdir:
            pathlib.Path(workdir).mkdir(parents=True, exist_ok=True)
            self._base_path = str(workdir)
        else:
            self._base_path = self.original_path or os.getcwd()

    def resolved_input(self):
        """A copy of ``input`` with the lattice and the files it names as absolute paths,
        resolved from the directory the input came from."""
        if self._input is None:
            raise ValueError("there is no input: give an Input, a path or TOML text")
        inp = _copy.deepcopy(self._input)
        base = pathlib.Path(self.original_path or os.getcwd())
        for key in PATH_KEYS:
            v = inp[key]
            if v and not pathlib.Path(v).is_absolute():
                inp[key] = str((base / v).resolve())
        inp.field.file = [str((base / f).resolve()) if f and not pathlib.Path(f).is_absolute() else f
                          for f in inp.field.file]
        return inp

    def write_input(self, input_filename=None):
        """Write the resolved input file into the run's directory, or to
        ``input_filename``, and return its path."""
        path = pathlib.Path(input_filename) if input_filename else pathlib.Path(self._base_path) / self._input_name()
        return self.resolved_input().write(path)

    def configure(self):
        """Make the run's directory and write the input file there, with the initial
        particles as an openPMD-beamphysics file when there are any."""
        self.setup_workdir(self._workdir)
        if self._initial_particles is not None:
            beam = pathlib.Path(self._base_path) / "initial_particles.h5"
            self._initial_particles.write(str(beam))
            self._input.beam.position_file = str(beam)
            self._input.beam.file = ""
        self.write_input()
        self._configured = True

    def get_run_script(self, write_to_path=True):
        """The command that runs the input file, as a list, run from ``path``. With
        ``write_to_path`` it is also written there as a file named run."""
        from lucifer.runner import executable
        script = [str(executable(self._command)), self._input_name()]
        if self.threads is not None:
            script = [f"OMP_NUM_THREADS={int(self.threads)}"] + script
        if write_to_path and self._base_path:
            (pathlib.Path(self._base_path) / "run").write_text(" ".join(script) + "\n")
        return script

    def run(self, raise_on_error=True):
        """Configure if needed, write the input file as it now stands, run the program, and
        return the output, a Result. A refusal or a failed run raises LuciferError with the
        program's message, or with ``raise_on_error`` off sets ``error`` and returns the
        output the run left."""
        from lucifer.runner import LuciferError, Result, run
        # The input is written on every run, so a key changed since the last run is the
        # one this run reads.
        if not self._configured:
            self.configure()
        else:
            self.write_input()
        # The previous run's stats file may be open for reading, and the program cannot
        # write over a file a reader holds open.
        if self._output is not None:
            self._output.close()
        self._finished = self._error = False
        try:
            self._output = run(self._input_name(), cwd=self._base_path, exe=self._command,
                               threads=self.threads, timeout=self._timeout)
        except LuciferError as err:
            self._error = True
            self._output = Result(err.returncode, err.output, pathlib.Path(self.input_file),
                                  pathlib.Path(self._base_path), self._input.run.out_root)
            if raise_on_error:
                raise
        finally:
            self._finished = True
        self.vprint(self._output.output)
        return self._output

    def load_output(self, **kwargs):
        """Read the output a run left in ``path``, a run made outside this object included,
        and return it."""
        from lucifer.runner import Result
        out_root = self._input.run.out_root
        stats = pathlib.Path(self._base_path) / f"{out_root}.stats.h5"
        if not stats.is_file():
            raise FileNotFoundError(f"no stats file at {stats}")
        self._output = Result(0, "", pathlib.Path(self.input_file), pathlib.Path(self._base_path), out_root)
        self._finished = True
        return self._output

    def vprint(self, *args, **kwargs):
        """Print only with ``verbose``."""
        if self._verbose:
            print(*args, **kwargs)

    # ------------------------------------------------------------------
    # Identity, copying, archiving.

    def fingerprint(self):
        """A hash of the input file with every key written out, which two objects share
        exactly when they would run the same written input. An object restored from an
        archive holds the input as the run read it, with every value the program derived
        written out, so its fingerprint differs from one that left those to the program."""
        text = self.resolved_input().toml(explicit=True) if self._input is not None else ""
        return hashlib.sha256(text.encode()).hexdigest()

    def copy(self):
        """A copy with its own input and output, not configured when it ran in a temporary
        directory."""
        other = _copy.copy(self)
        other._input = _copy.deepcopy(self._input)
        # A copy of the output without its open stats file, which cannot be copied.
        other._output = dataclasses.replace(self._output) if self._output is not None else None
        other._initial_particles = _copy.deepcopy(self._initial_particles)
        other._tempdir = None
        other.reset()
        return other

    def reset(self):
        """Forget the run's directory when it was temporary, so the next run configures a
        new one."""
        if self._use_temp_dir:
            self._base_path = None
            self._configured = False

    def archive(self, h5=None):
        """Copy the run's stats file, which holds the input file in meta/input_echo and every
        output array, to ``h5``, by default lucifer_<fingerprint>.h5 here, and return its
        path."""
        if self._output is None or not self._output.stats_path.is_file():
            raise ValueError("there is nothing to archive: run() first")
        path = pathlib.Path(h5 or f"lucifer_{self.fingerprint()[:12]}.h5")
        shutil.copyfile(self._output.stats_path, path)
        return path

    def load_archive(self, h5, configure=True):
        """Take the input file from an archive's meta/input_echo and its output from the
        archive itself."""
        import h5py
        from lucifer.models import Input
        from lucifer.runner import Result
        h5 = pathlib.Path(h5).resolve()
        with h5py.File(h5, "r") as f:
            echo = f["meta/input_echo"][()]
        self._input = Input.loads(echo.decode() if isinstance(echo, bytes) else str(echo))
        self.original_path = str(h5.parent)
        self._archive = h5.name
        self._output = Result(0, "", h5, h5.parent, self._input.run.out_root, stats_file=h5)
        self._finished = True
        self._configured = False

    @classmethod
    def from_archive(cls, archive_h5):
        """A Lucifer rebuilt from an archive: its input from meta/input_echo, its output the
        archive."""
        c = cls()
        c.load_archive(archive_h5)
        return c

    def to_hdf5(self, filename):
        """archive() to ``filename``."""
        return self.archive(filename)

    @classmethod
    def from_hdf5(cls, filename):
        """from_archive() of ``filename``."""
        return cls.from_archive(filename)

    # ------------------------------------------------------------------
    # Plotting and showing.

    def _need_output(self):
        if self._output is None:
            raise ValueError("there is no output: run() first")
        return self._output

    def plot(self, y="field/total/power", **kwargs):
        """lucifer.plot.plot of the output's stats file, which names its arguments."""
        from lucifer.plot import plot
        return plot(self._need_output().stats_path, y=y, **kwargs)

    def summary(self, **kwargs):
        """lucifer.plot.summary of the output's stats file."""
        from lucifer.plot import summary
        return summary(self._need_output().stats_path, **kwargs)

    def _state(self):
        if self._error:
            return "error"
        if self._finished:
            return "finished"
        return "configured" if self._configured else "not run"

    def __repr__(self):
        if self._archive:
            return f"Lucifer(archive={self._archive!r}, {self._state()})"
        name = self.original_input_file or ("in memory" if self._input is not None else None)
        return f"Lucifer(input={name!r}, path={self._base_path!r}, {self._state()})"

    def _repr_html_(self):
        from lucifer._display import short, table
        rows = [("archive", self._archive)] if self._archive else []
        rows += [("input file", self.input_file or "in memory"), ("path", self._base_path or "not configured"),
                 ("state", self._state())]
        if self._input is not None:
            rows += [(k, short(v)) for k, v in self._input.items()]
        if self._output is not None:
            rows += [("stats file", self._output.stats_path.name)]
        return table("Lucifer", rows)
