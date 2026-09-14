"""A small process pool for a check script's independent scan points.

The check scripts run their scan points as separate tracker (or genesis4)
processes, historically one after another. The points are independent, so
running a few at once changes NO physics, NO configuration and NO tolerance --
only wall time. Each job callable launches one process and blocks on it;
threads suffice because the work happens in the child process.
"""

from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor


def run_all(jobs, threads_per_job=4):
    """Run the job callables a few at a time and return their results in order.

    The worker count is the caller's core allowance over threads_per_job, the OMP
    thread count each child runs with. The allowance is LUCIFER_CPU_BUDGET when the
    caller set it, which the benchmark does for every check section from its own
    share, so several sections at once cannot each take the whole machine. Unset,
    which is a script run by hand, the allowance is the machine. The first failure
    propagates (a job that calls sys.exit raises SystemExit out of its future) after
    all jobs were launched.
    """
    jobs = list(jobs)
    if not jobs:
        return []
    budget = os.environ.get("LUCIFER_CPU_BUDGET", "").strip()
    cpus = int(budget) if budget else (os.cpu_count() or 8)
    workers = max(1, cpus // max(1, int(threads_per_job)))
    workers = min(workers, len(jobs))
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = [ex.submit(job) for job in jobs]
        return [f.result() for f in futures]
