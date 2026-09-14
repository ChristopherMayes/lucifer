"""Hold a lock file while a command runs.

    with_lock.py <lockfile> <command> [args...]

The lock is an exclusive flock on the named file, taken before the command starts and
held until it ends. The kernel releases it when the process dies, however it dies, so a
killed run leaves no stale lock behind, which a lock directory would. The command's
stdin, stdout and stderr are inherited, and its exit status is returned.

The benchmark's device section runs under this. The two build passes of a keystone run
at once, and the device is one resource across both, so their device sections take
turns while their CPU sections carry on (doc/validation.md). The wait, when there was
one, is printed so the cost of taking turns is a measured number.
"""

import fcntl
import os
import subprocess
import sys
import time


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    lockfile, cmd = sys.argv[1], sys.argv[2:]
    os.makedirs(os.path.dirname(os.path.abspath(lockfile)), exist_ok=True)
    t0 = time.monotonic()
    with open(lockfile, "a") as lf:
        fcntl.flock(lf, fcntl.LOCK_EX)
        waited = time.monotonic() - t0
        if waited > 1.0:
            print(f"  device lock: waited {waited:.0f} s for {lockfile}", flush=True)
        else:
            print(f"  device lock: taken at once, {lockfile}", flush=True)
        r = subprocess.run(cmd)
    sys.exit(r.returncode)


if __name__ == "__main__":
    main()
