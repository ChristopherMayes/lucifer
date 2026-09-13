"""The keystone as one pytest suite (doc/validation.md, val-the-keystone-rule).

Five jobs run once, at the same time, each from an absolute working directory with its
own log, and every outcome is kept: a zero return, a nonzero return, a launch that
failed, a job that timed out. Five tests read those outcomes, one each. A sixth test
regenerates the three page sets and requires the Markdown diff to be empty, and it
verifies for itself that every one of the five outcomes is present and successful
before it runs a generator, so it does not lean on test order and cannot pass on a
partial run. Nothing here builds: both trees are prerequisites and are checked before
any job launches, since a benchmark reads the source tree throughout.

    BUILD_PRODUCTION=N ./util/conda_compile
    ./util/conda_compile
    <bmad-fel-validate python> -m pytest lucifer/tests/test_keystone.py -v

What the contract means in practice:

- Selecting one job test still runs all five jobs. The fixture is session scoped and the
  keystone is whole or it is nothing, so the cost of the smallest selection is the cost
  of the keystone.
- xdist is refused (conftest.py beside this file). A session fixture runs once per
  worker, so a parallel run would launch the five jobs several times over, into the
  same work directories.
- The regeneration test writes into the tree: the three generated page sets under
  lucifer/doc/generated. That is the keystone's own contract, and a clean tree is what
  an empty diff then proves.
- A job that runs past its timeout is killed as a process group. Each job is started in
  its own session, so the group is the shell wrapper and everything it launched, the
  tracker binaries and the check scripts included. A timeout on the wrapper alone would
  leave those running against the work directory of the next keystone.
- The artifacts of a session, the five logs, the two results files and the two work
  directories, live under one fresh directory named in the session header, or under
  KEYSTONE_ARTIFACTS when that is set. They are kept, since a failure is read from them.

The five check scripts are unchanged. This file launches them the way the recorded
recipe did and replaces the retyping of that recipe, which twice launched a benchmark
from a directory an earlier command had moved into, ran four jobs of five, and
reported a complete keystone from a one-line log.
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

ROOT = pathlib.Path(__file__).resolve().parents[2]
DEBUG_EXE = ROOT / "debug" / "bin" / "lucifer"
PRODUCTION_EXE = ROOT / "production" / "bin" / "lucifer"
BENCHMARK = ROOT / "lucifer" / "tests" / "run_fel_benchmark.sh"
WAVEFRONT = ROOT / "lucifer" / "wavefront" / "tests" / "run_validation.sh"
EXAMPLES = ROOT / "lucifer" / "examples" / "run_examples.sh"
SCRIPTS = ROOT / "lucifer" / "tests" / "scripts"
GENERATED = ROOT / "lucifer" / "doc" / "generated"

# One hour a job. The keystone takes nine minutes cold, so a job still running at an
# hour has hung, and the group is killed rather than waited on.
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


def jobs(art: pathlib.Path) -> tuple[Job, ...]:
    """The five jobs as the recorded recipe launches them, with this session's artifacts."""
    return (
        Job("debug", (str(BENCHMARK), "--results", str(art / "fel-debug.txt"),
                      "--work-dir", str(art / "wd-dbg")), ROOT),
        Job("production", (str(BENCHMARK), "--exe", str(PRODUCTION_EXE),
                           "--results", str(art / "fel-prod.txt"),
                           "--work-dir", str(art / "wd-prd")), ROOT),
        Job("regression", (sys.executable, "-m", "pytest", "test_fortran.py",
                           f"--bmad-bin={ROOT / 'debug' / 'bin'}"), ROOT / "regression_tests"),
        Job("wavefront", (str(WAVEFRONT),), ROOT),
        # One worker for the examples, which is this job's whole allocation while the
        # five run together. The examples are not on the critical path, and every worker
        # they take is taken from the two benchmark passes that are: measured on
        # 2026-09-13, the examples at three workers finished in 347 s against 566 s and
        # 570 s at one, and the debug pass they share the machine with went from 692 s
        # and 694 s to 748 s.
        # A standalone run of run_examples.sh takes its own default instead.
        Job("examples", (str(EXAMPLES), "--no-figures", "--jobs", "1"), ROOT),
    )


def run_job(job: Job, log: pathlib.Path) -> Outcome:
    """Run one job to its return code, or to a recorded failure, never to an exception."""
    out = Outcome(job.name, log)
    t0 = time.monotonic()
    env = dict(os.environ, **dict(job.env))
    try:
        with open(log, "w") as lf:
            proc = subprocess.Popen(list(job.argv), cwd=str(job.cwd), env=env,
                                    stdout=lf, stderr=subprocess.STDOUT,
                                    start_new_session=True)
            try:
                out.returncode = proc.wait(timeout=job.timeout)
            except subprocess.TimeoutExpired:
                # The wrapper's whole group, not the wrapper alone.
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
                    + ", ".join(missing))
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
    """Launch the five jobs at once and keep every outcome. Not a success gate."""
    preflight()
    chosen = os.environ.get("KEYSTONE_ARTIFACTS")
    art = pathlib.Path(chosen).resolve() if chosen else \
        pathlib.Path(tempfile.mkdtemp(prefix="lucifer-keystone-"))
    art.mkdir(parents=True, exist_ok=True)
    print(f"\nkeystone: five jobs from {ROOT}, artifacts in {art}")
    outcomes: dict[str, Outcome] = {}
    with ThreadPoolExecutor(max_workers=len(JOB_NAMES)) as pool:
        futures = {j.name: pool.submit(run_job, j, art / f"{j.name}.log") for j in jobs(art)}
        for name, fut in futures.items():
            try:
                outcomes[name] = fut.result()
            except Exception as exc:                   # run_job does not raise; belt and braces
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
    """Regenerate the three page sets and require no diff, after all five jobs passed.

    The prerequisite is checked here and not inherited from the five tests above, since
    those may be deselected or ordered after this one. A partial or failed run fails
    this test; it never skips it.
    """
    bad = [n for n in JOB_NAMES if not keystone.outcomes.get(n, Outcome(n, pathlib.Path())).passed]
    if bad:
        pytest.fail("regeneration refused: these jobs did not pass: " + ", ".join(bad)
                    + f" (artifacts in {keystone.artifacts})")
    art = keystone.artifacts
    generators = (
        (str(SCRIPTS / "report_validation.py"), "--debug", str(art / "fel-debug.txt"),
         "--production", str(art / "fel-prod.txt"), "--out", str(GENERATED / "validation-measured.md")),
        (str(SCRIPTS / "report_api.py"), "--code", str(ROOT / "lucifer" / "code"),
         "--code", str(ROOT / "lucifer" / "program"), "--out", str(GENERATED / "api.md")),
        (str(SCRIPTS / "report_examples.py"), "--examples", str(ROOT / "lucifer" / "examples"),
         "--out", str(GENERATED / "examples")),
    )
    for gen in generators:
        r = subprocess.run([sys.executable, *gen], cwd=str(ROOT), capture_output=True, text=True)
        assert r.returncode == 0, f"{pathlib.Path(gen[0]).name} exited {r.returncode}\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}"
    # The Markdown only. The figures are committed and a plotting library upgrade
    # rewrites them all with no physics moving (doc/validation.md).
    r = subprocess.run(["git", "diff", "--exit-code", "--stat", "--",
                        "lucifer/doc/generated/*.md", "lucifer/doc/generated/examples/*.md"],
                       cwd=str(ROOT), capture_output=True, text=True)
    assert r.returncode == 0, ("a generated page moved after regeneration; a moved digit is a bug, "
                               "not a new baseline:\n" + r.stdout)
