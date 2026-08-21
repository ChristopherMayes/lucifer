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

    The worker count is sized so the machine is fully but not oversubscribed:
    cpu_count / threads_per_job, where threads_per_job is the OMP thread count
    each child runs with. The first failure propagates (a job that calls
    sys.exit raises SystemExit out of its future) after all jobs were launched.
    """
    jobs = list(jobs)
    if not jobs:
        return []
    workers = max(1, (os.cpu_count() or 8) // max(1, int(threads_per_job)))
    workers = min(workers, len(jobs))
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = [ex.submit(job) for job in jobs]
        return [f.result() for f in futures]
