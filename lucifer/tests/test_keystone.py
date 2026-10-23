"""The keystone as one pytest suite (doc/validation.md, val-the-keystone-rule).

Five jobs run once, at the same time. Each runs from an absolute working directory with
its own log. Every outcome is kept: a zero return, a nonzero return, a launch that
failed, a job that timed out. Five tests read those outcomes, one each. A sixth test
regenerates the five page sets and requires the Markdown diff to be empty. Before it
runs a generator, the sixth test verifies for itself that all five outcomes are present
and successful. So it does not lean on test order and cannot pass on a partial run.
Nothing here builds. Both trees are prerequisites and are checked before any job
launches, since a benchmark reads the source tree throughout.

    BUILD_PRODUCTION=N lucifer/util/conda_compile
    lucifer/util/conda_compile
    <bmad-fel-validate python> -m pytest lucifer/tests/test_keystone.py -v

What the contract means in practice:

- Selecting one job test still runs all five jobs. The fixture is session scoped, so
  the smallest selection costs a whole keystone.
- xdist is refused (conftest.py beside this file). A session fixture runs once per
  worker, so a parallel run would launch the five jobs several times over, into the
  same work directories.
- The regeneration test writes into the tree: the four generated page sets under
  lucifer/doc/generated. The keystone writes them by design. An empty diff then shows
  that the committed pages agree with what the sources generate.
- A job that runs past its timeout is killed as a process group. Each job is started in
  its own session, so the group is the shell wrapper and everything it launched, the
  tracker binaries and the check scripts included. A timeout on the wrapper alone would
  leave those running against the work directory of the next keystone.
- The artifacts of a session (the five logs, the two results files and the two work
  directories) live under one fresh directory named in the session header. Where
  KEYSTONE_ARTIFACTS is set they live there instead. They are kept, since a failure is
  read from them.

The five check scripts are unchanged. This file launches them as the recipe in
doc/validation.md does. A retyped recipe twice launched a benchmark from a directory an
earlier command had moved into, ran four jobs of five, and reported a complete keystone
from a one-line log. This file replaces the retyping.
"""

from __future__ import annotations

import dataclasses
import os
import pathlib
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor

import pytest

# Four locations. Each defaults to the layout of a Lucifer source tree inside the Bmad
# checkout, which then runs with nothing set. This is the rule of util/locations.sh,
# which the shell scripts source.
#   The Lucifer source is this file's tree.
#   The Bmad source, with the regression suite and its debug binaries, is BMAD_DIR,
#     else DIST_BASE_DIR as util/dist_source_me exports it, else the Lucifer source's parent.
#   The two binary directories default to debug/ and production/ beside the Lucifer
#     source, where Master.cmake places a project's binaries. That holds in the Bmad
#     checkout and out of it. LUCIFER_DEBUG_BIN and LUCIFER_PRODUCTION_BIN name them elsewhere.
# A relative value is taken against the directory pytest was started in. Every job is
# handed the absolute values, since the jobs run from the Lucifer source and a relative
# path would mean something else there.
LUCIFER = pathlib.Path(__file__).resolve().parents[1]


def bmad_source() -> tuple[pathlib.Path, list[str]]:
    """The Bmad checkout, and the notes on any variable passed over.

    A variable naming a directory with no util/dist_source_me is passed over. The
    reason: dist_source_me uses BMAD_DIR as a loop variable and leaves it at the
    checkout's util/ once sourced. DIST_BASE_DIR, which it exports, still names the
    checkout.
    """
    notes = []
    for var in ("BMAD_DIR", "DIST_BASE_DIR"):
        value = os.environ.get(var)
        if not value:
            continue
        path = pathlib.Path(os.path.abspath(value))
        if (path / "util" / "dist_source_me").is_file():
            return path, notes
        notes.append(f"{var}={value} holds no util/dist_source_me and is passed over")
    return LUCIFER.parent, notes


BMAD, BMAD_NOTES = bmad_source()
DEBUG_BIN = pathlib.Path(os.path.abspath(os.environ.get("LUCIFER_DEBUG_BIN") or LUCIFER.parent / "debug" / "bin"))
PRODUCTION_BIN = pathlib.Path(os.path.abspath(os.environ.get("LUCIFER_PRODUCTION_BIN") or LUCIFER.parent / "production" / "bin"))
# What every job is handed: each value as it was given, made absolute, with no symbolic
# link followed. openPMD-beamphysics is not a location. It is a package installed in the
# interpreter's environment, which every job inherits.
LOCATIONS = {"BMAD_DIR": str(BMAD),
             "LUCIFER_DEBUG_BIN": str(DEBUG_BIN), "LUCIFER_PRODUCTION_BIN": str(PRODUCTION_BIN)}
DEBUG_EXE = DEBUG_BIN / "lucifer"
PRODUCTION_EXE = PRODUCTION_BIN / "lucifer"
BENCHMARK = LUCIFER / "tests" / "run_fel_benchmark.sh"
WAVEFRONT = LUCIFER / "wavefront" / "tests" / "run_validation.sh"
EXAMPLES = LUCIFER / "examples" / "run_examples.sh"
SCRIPTS = LUCIFER / "tests" / "scripts"
GENERATED = LUCIFER / "doc" / "generated"

# One hour a job. The keystone takes nine minutes cold, so a job still running at an
# hour has hung. Its process group is then killed rather than waited on.
TIMEOUT_S = 3600

JOB_NAMES = ("debug", "production", "regression", "wavefront", "examples")


@dataclasses.dataclass(frozen=True)
class Job:
    """One keystone job: what to run, where, with what, and for how long at most."""
    name: str
    argv: tuple[str, ...]
    cwd: pathlib.Path
    env: tuple[tuple[str, str], ...] = ()
    timeout: float = TIMEOUT_S


@dataclasses.dataclass
class Outcome:
    """What one job did. error is set when the job did not run to a return code."""
    name: str
    log: pathlib.Path
    returncode: int | None = None
    error: str | None = None
    seconds: float = 0.0

    @property
    def passed(self) -> bool:
        return self.error is None and self.returncode == 0


@dataclasses.dataclass
class Keystone:
    artifacts: pathlib.Path
    outcomes: dict[str, Outcome]


def performance_cores() -> int:
    """The cores the static allocation in jobs() divides. These are the performance
    cores on a Mac and the logical cores elsewhere."""
    for args in (["sysctl", "-n", "hw.perflevel0.physicalcpu"], ["sysctl", "-n", "hw.ncpu"]):
        try:
            n = int(subprocess.run(args, capture_output=True, text=True, check=True).stdout)
            if n > 0:
                return n
        except (OSError, ValueError, subprocess.CalledProcessError):
            pass
    return os.cpu_count() or 4


def jobs(art: pathlib.Path) -> tuple[Job, ...]:
    """The five jobs as the recipe in doc/validation.md launches them, with this session's artifacts.

    The static allocation. The machine's performance cores are divided once, here, and
    every job spends what it is handed and no more.
      Each benchmark pass is handed the whole count in --cpus. Its check sections run
        cores/4 at once, and each section's pools take the cores over that number
        (run_fel_benchmark.sh).
      The two passes' device sections share one lock, since the device is one resource.
      The examples take two workers.
      The regression suite and the wavefront validation take one process each.
    The sum at the widest is two passes at the full count, two example runs at the
    thread count their runtime chooses, and two more processes. That sum is above the
    core count by choice. doc/validation.md (val-the-keystone-rule) records the wall
    times that set it.
    """
    cores = performance_cores()
    lock = (("LUCIFER_DEVICE_LOCK", str(art / "device.lock")),)
    # The two passes reach their tier blocks together. Eleven single-thread processes a
    # pass make twenty-two on twelve cores. The debug pass, which finishes last, runs
    # its tiers first and then writes a marker. The production pass holds its own tier
    # block until that marker appears and runs the rest of its work meanwhile.
    tiers = str(art / "tiers-debug.done")
    return (
        Job("debug", (str(BENCHMARK), "--exe", str(DEBUG_EXE), "--build", "debug",
                      "--results", str(art / "fel-debug.txt"),
                      "--work-dir", str(art / "wd-dbg"), "--cpus", str(cores),
                      "--tiers-marker", tiers), LUCIFER, lock),
        Job("production", (str(BENCHMARK), "--exe", str(PRODUCTION_EXE), "--build", "production",
                           "--results", str(art / "fel-prod.txt"),
                           "--work-dir", str(art / "wd-prd"), "--cpus", str(cores),
                           "--wait-tiers", tiers), LUCIFER, lock),
        Job("regression", (sys.executable, "-m", "pytest", "test_fortran.py",
                           f"--bmad-bin={BMAD / 'debug' / 'bin'}"), BMAD / "regression_tests"),
        Job("wavefront", (str(WAVEFRONT),), LUCIFER),
        # Two workers for the examples, this job's whole allocation while the five run
        # together. At one worker the examples job finishes last. At two the debug pass
        # finishes last. Three workers slowed the benchmark passes by more than they
        # saved. doc/validation.md (val-the-keystone-rule) records the wall times. A
        # standalone run of run_examples.sh takes its own default instead.
        Job("examples", (str(EXAMPLES), "--no-figures", "--jobs", "2"), LUCIFER),
    )


def run_job(job: Job, log: pathlib.Path) -> Outcome:
    """Run one job to its return code, or to a recorded failure, never to an exception."""
    out = Outcome(job.name, log)
    t0 = time.monotonic()
    env = dict(os.environ, **LOCATIONS, **dict(job.env))
    try:
        with open(log, "w") as lf:
            proc = subprocess.Popen(list(job.argv), cwd=str(job.cwd), env=env,
                                    stdout=lf, stderr=subprocess.STDOUT,
                                    start_new_session=True)
            try:
                out.returncode = proc.wait(timeout=job.timeout)
            except subprocess.TimeoutExpired:
                # Kill the wrapper's whole process group, and not the wrapper alone.
                try:
                    os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                except ProcessLookupError:
                    pass
                proc.wait()
                out.error = f"timed out after {job.timeout:.0f} s; process group killed"
    except Exception as exc:                       # a launch that never became a process
        out.error = f"launch failed: {type(exc).__name__}: {exc}"
    out.seconds = time.monotonic() - t0
    return out


def tail(log: pathlib.Path, n: int = 40) -> str:
    try:
        lines = log.read_text(errors="replace").splitlines()
    except OSError as exc:
        return f"(log unreadable: {exc})"
    return "\n".join(lines[-n:])


def preflight() -> None:
    """What must be true before a single job launches."""
    missing = [str(p) for p in (DEBUG_EXE, PRODUCTION_EXE) if not os.access(p, os.X_OK)]
    if missing:
        pytest.fail("both builds are prerequisites and the suite never builds; missing: "
                    + ", ".join(missing) + " (LUCIFER_DEBUG_BIN and LUCIFER_PRODUCTION_BIN name "
                    "the binary directories where they are not beside the Lucifer source)")
    if not (BMAD / "regression_tests" / "test_fortran.py").is_file() or not (BMAD / "debug" / "bin").is_dir():
        pytest.fail(f"no Bmad source with its regression suite and a debug build at {BMAD} "
                    "(BMAD_DIR names the Bmad checkout where it is not the Lucifer source's parent)")
    # The release comparison's own cases, then the Bmad found against lucifer/BMAD_RELEASE,
    # the same routine the build script and the Linux CI job call.
    check = LUCIFER / "util" / "bmad_release_check"
    for arg in ("--self-test", str(BMAD)):
        r = subprocess.run([str(check), arg], capture_output=True, text=True)
        for line in r.stdout.splitlines():
            print(f"keystone: {line}")
        if r.returncode != 0:
            pytest.fail(f"bmad_release_check {arg} refused:\n{r.stdout}{r.stderr}")
    missing = [str(p) for p in (BENCHMARK, WAVEFRONT, EXAMPLES) if not os.access(p, os.X_OK)]
    if missing:
        pytest.fail("keystone scripts not executable: " + ", ".join(missing))
    if "bmad-fel-validate" not in sys.executable:
        pytest.fail("run this suite with the bmad-fel-validate interpreter, which runs the "
                    f"regression suite and the generators; got {sys.executable}")
    if shutil.which("git") is None:
        pytest.fail("git is needed for the regeneration diff")


@pytest.fixture(scope="session")
def keystone() -> Keystone:
    """Launch the five jobs at once and keep every outcome.

    A failed job does not fail the fixture. The tests read the outcomes. A failed
    preflight does fail it, before any job launches.
    """
    preflight()
    chosen = os.environ.get("KEYSTONE_ARTIFACTS")
    art = pathlib.Path(chosen).resolve() if chosen else \
        pathlib.Path(tempfile.mkdtemp(prefix="lucifer-keystone-"))
    art.mkdir(parents=True, exist_ok=True)
    print(f"\nkeystone: five jobs, Lucifer {LUCIFER}, Bmad {BMAD}, artifacts in {art}")
    for note in BMAD_NOTES:
        print(f"keystone: {note}")
    outcomes: dict[str, Outcome] = {}
    with ThreadPoolExecutor(max_workers=len(JOB_NAMES)) as pool:
        futures = {j.name: pool.submit(run_job, j, art / f"{j.name}.log") for j in jobs(art)}
        for name, fut in futures.items():
            try:
                outcomes[name] = fut.result()
            except Exception as exc:                   # run_job does not raise. A pool fault lands here.
                outcomes[name] = Outcome(name, art / f"{name}.log",
                                         error=f"runner failed: {type(exc).__name__}: {exc}")
    for name in JOB_NAMES:
        o = outcomes[name]
        state = "passed" if o.passed else (o.error or f"exit {o.returncode}")
        print(f"keystone: {name:11s} {o.seconds:7.1f} s  {state}")
    return Keystone(art, outcomes)


@pytest.mark.parametrize("name", JOB_NAMES)
def test_job(keystone: Keystone, name: str) -> None:
    o = keystone.outcomes.get(name)
    assert o is not None, f"{name}: no outcome recorded"
    assert o.passed, (f"{name}: {o.error or f'exit {o.returncode}'} after {o.seconds:.0f} s\n"
                      f"log: {o.log}\n{tail(o.log)}")


def test_regeneration(keystone: Keystone) -> None:
    """Regenerate the five page sets and require no diff, after all five jobs passed.

    The prerequisite is checked here and not inherited from the five tests above, since
    those may be deselected or ordered after this one. A partial or failed run fails
    this test. It never skips it.
    """
    bad = [n for n in JOB_NAMES if not keystone.outcomes.get(n, Outcome(n, pathlib.Path())).passed]
    if bad:
        pytest.fail("regeneration refused: these jobs did not pass: " + ", ".join(bad)
                    + f" (artifacts in {keystone.artifacts})")
    art = keystone.artifacts
    generators = (
        (str(SCRIPTS / "report_validation.py"), "--debug", str(art / "fel-debug.txt"),
         "--production", str(art / "fel-prod.txt"), "--out", str(GENERATED / "validation-measured.md")),
        (str(SCRIPTS / "report_api.py"), "--code", str(LUCIFER / "code"),
         "--code", str(LUCIFER / "program"), "--out", str(GENERATED / "api.md")),
        (str(SCRIPTS / "report_examples.py"), "--examples", str(LUCIFER / "examples"),
         "--out", str(GENERATED / "examples")),
        # The startup-noise page's tables, from the committed JSON of the driver's runs. The
        # driver runs on demand, so this costs no run.
        (str(SCRIPTS / "report_startup_noise.py"), "--json", str(GENERATED / "startup-noise"),
         "--out", str(GENERATED / "startup-noise-measured.md")),
        # The layouts page runs its own fixtures with the production executable. A failed
        # fixture fails this test, keeps its log in the artifact directory and leaves the
        # committed page as it was, since the generator writes the page whole or not at all.
        (str(SCRIPTS / "report_layouts.py"), "--exe", str(PRODUCTION_EXE),
         "--out", str(GENERATED / "layouts.md"), "--log", str(art / "layouts.log")),
    )
    for gen in generators:
        r = subprocess.run([sys.executable, *gen], cwd=str(LUCIFER), capture_output=True, text=True)
        assert r.returncode == 0, f"{pathlib.Path(gen[0]).name} exited {r.returncode}\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}"
    # The Markdown only. The figures are committed and a plotting library upgrade
    # rewrites them all with no physics moving (doc/validation.md).
    r = subprocess.run(["git", "diff", "--exit-code", "--stat", "--",
                        "doc/generated/*.md", "doc/generated/examples/*.md", "examples/measured.tsv"],
                       cwd=str(LUCIFER), capture_output=True, text=True)
    assert r.returncode == 0, ("a generated page moved after regeneration; a moved digit is a bug, "
                               "not a new baseline:\n" + r.stdout)
