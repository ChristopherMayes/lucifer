#!/usr/bin/env python3
"""Measure the FFTW plans of the transverse transform and write them as a Fortran source.

    util/make_fft_wisdom.py [--lib <libfftw3>] [--sessions 3] [--out src/wavefront_wisdom_mod.f90]

wavefront_fft2 (src/wavefront_mod.f90) transforms one n by n plane in place, forward and
backward, on a buffer of fftw_alloc_complex. FFTW's planner finds a fast plan for that
problem by timing candidates, and the plan it finds differs from one planning to the
next. This script plans once, at FFTW_PATIENT, for every power-of-two n from 32 to 2048,
and writes FFTW's record of the plans it chose, its wisdom, into the module
wavefront_wisdom_mod. The planner of wavefront_fft2 then makes those same plans in every
run, from the text and with no timing.

One session plans every size in turn, the wisdom accumulating, so the plans of one
session agree with each other where they share a part. Several sessions are run and the
one whose transforms are fastest is kept. The kept text is then read back into an empty
planner, and each size is planned from it alone and timed beside the estimated plan.

The wisdom belongs to the FFTW version and build it was measured with. Another version
does not accept it, and wavefront_fft2 then plans with FFTW_ESTIMATE. Run this script on
a machine doing nothing else, with the library the build links, and commit the module.
Every recorded digit of the validation suite moves with it, in the last places.
"""
import argparse
import ctypes
import ctypes.util
import datetime
import os
import platform
import sys
import time

SIZES = [32, 64, 128, 256, 512, 1024, 2048]
FFTW_FORWARD, FFTW_BACKWARD = -1, 1
FFTW_ESTIMATE, FFTW_PATIENT, FFTW_WISDOM_ONLY = 1 << 6, 1 << 5, 1 << 21
LINE = 120     # The longest line of wisdom text in one Fortran string.


def find_library():
    """Return the path of libfftw3 in the active conda environment, or what the loader finds."""
    prefix = os.environ.get("CONDA_PREFIX", "")
    for name in ("libfftw3.dylib", "libfftw3.so"):
        path = os.path.join(prefix, "lib", name)
        if prefix and os.path.exists(path):
            return path
    return ctypes.util.find_library("fftw3")


def bind(path):
    """Load libfftw3 and declare the entry points this script calls."""
    lib = ctypes.CDLL(path)
    vp, ci, cu = ctypes.c_void_p, ctypes.c_int, ctypes.c_uint
    lib.fftw_alloc_complex.restype, lib.fftw_alloc_complex.argtypes = vp, [ctypes.c_size_t]
    lib.fftw_free.restype, lib.fftw_free.argtypes = None, [vp]
    lib.fftw_plan_dft_2d.restype, lib.fftw_plan_dft_2d.argtypes = vp, [ci, ci, vp, vp, ci, cu]
    lib.fftw_execute.restype, lib.fftw_execute.argtypes = None, [vp]
    lib.fftw_destroy_plan.restype, lib.fftw_destroy_plan.argtypes = None, [vp]
    lib.fftw_forget_wisdom.restype, lib.fftw_forget_wisdom.argtypes = None, []
    lib.fftw_export_wisdom_to_string.restype, lib.fftw_export_wisdom_to_string.argtypes = vp, []
    lib.fftw_import_wisdom_from_string.restype = ci
    lib.fftw_import_wisdom_from_string.argtypes = [ctypes.c_char_p]
    return lib


def plan_pair(lib, n, buf, flags):
    """Return the forward and backward in-place plans of an n by n transform, or None."""
    fwd = lib.fftw_plan_dft_2d(n, n, buf, buf, FFTW_FORWARD, flags)
    bwd = lib.fftw_plan_dft_2d(n, n, buf, buf, FFTW_BACKWARD, flags)
    if not fwd or not bwd:
        for plan in (fwd, bwd):
            if plan:
                lib.fftw_destroy_plan(plan)
        return None
    return fwd, bwd


def time_pair(lib, n, buf, pair):
    """Return the seconds of one transform, the best of five timings of a pair."""
    reps = max(2, (1 << 22) // (n * n))
    ctypes.memset(buf, 0, 16 * n * n)
    best = float("inf")
    for _ in range(5):
        t0 = time.perf_counter()
        for _ in range(reps):
            lib.fftw_execute(pair[0])
            lib.fftw_execute(pair[1])
        best = min(best, (time.perf_counter() - t0) / (2 * reps))
    return best


def session(lib, flags):
    """Plan every size in turn under flags and return the seconds of a transform at each."""
    seconds = {}
    for n in SIZES:
        buf = lib.fftw_alloc_complex(n * n)
        pair = plan_pair(lib, n, buf, flags)
        if pair is None:
            sys.exit(f"no plan for {n} by {n}")
        seconds[n] = time_pair(lib, n, buf, pair)
        for plan in pair:
            lib.fftw_destroy_plan(plan)
        lib.fftw_free(buf)
    return seconds


def exported(lib):
    """Return the planner's wisdom as text."""
    ptr = lib.fftw_export_wisdom_to_string()
    text = ctypes.string_at(ptr).decode("ascii")
    ctypes.CDLL(None).free(ctypes.c_void_p(ptr))
    return text


def fortran_source(text, version, made_on, timings):
    """Return the module wavefront_wisdom_mod holding the wisdom text."""
    pieces = []
    for line in text.splitlines():
        line = line.rstrip()
        while len(line) > LINE:
            cut = line.rfind(" ", 0, LINE)
            pieces.append(line[:cut])
            line = line[cut:]
        pieces.append(line)
    if any("'" in p for p in pieces):
        sys.exit("the wisdom text holds a quote, which the Fortran strings do not carry")
    rows = ", &\n".join(f"    '{p}'" for p in pieces)
    table = "\n".join(f"!   {n:5d}   {timings[n][0] * 1e6:10.1f}   {timings[n][1] * 1e6:10.1f}"
                      for n in SIZES)
    return f"""!+
! Module wavefront_wisdom_mod
!
! The FFTW wisdom of the transverse transform, written by util/make_fft_wisdom.py.
! Do not edit this file: run the script and commit what it writes.
!
! wavefront_fft2 (wavefront_mod) transforms an n by n plane in place. The planner of
! FFTW finds a fast plan for that by timing, and what it finds differs from one planning
! to the next. This text is its record of the plans one planning chose, for every
! power-of-two n from {SIZES[0]} to {SIZES[-1]}. wavefront_fft2_plan reads it and makes the same plans
! in every run. A size it does not hold, or an FFTW that does not accept it, is planned
! with FFTW_ESTIMATE.
!
! Made with {version}, {made_on}.
! One transform on one thread, in microseconds, as the script timed it:
!
!       n     estimated       wisdom
{table}
!-

module wavefront_wisdom_mod

implicit none

! The lines of the wisdom text. A line of FFTW's is one or more elements here, and
! wavefront_fft2_plan joins them.

character(*), parameter :: wf_wisdom_made_with = '{version}'
character({LINE}), parameter :: wf_wisdom_line(*) = [character({LINE}) :: &
{rows}]

end module
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--lib", default=find_library(), help="libfftw3 to measure with")
    ap.add_argument("--sessions", type=int, default=3, help="planning sessions, the fastest kept")
    ap.add_argument("--out", default="src/wavefront_wisdom_mod.f90", help="the module to write")
    args = ap.parse_args()
    if not args.lib:
        sys.exit("no libfftw3 found: give --lib")
    lib = bind(args.lib)
    version = (ctypes.c_char * 64).in_dll(lib, "fftw_version").value.decode("ascii")

    lib.fftw_forget_wisdom()
    estimated = session(lib, FFTW_ESTIMATE)

    kept = None
    for k in range(args.sessions):
        lib.fftw_forget_wisdom()
        seconds = session(lib, FFTW_PATIENT)
        total = sum(seconds[n] / estimated[n] for n in SIZES)
        print(f"session {k + 1}: " + "  ".join(f"{n}: {seconds[n] * 1e6:.1f}" for n in SIZES)
              + f"   sum of ratios {total:.3f}")
        if kept is None or total < kept[0]:
            kept = (total, exported(lib))

    # The kept text, read back into an empty planner and planned from alone.
    lib.fftw_forget_wisdom()
    if not lib.fftw_import_wisdom_from_string(kept[1].encode("ascii")):
        sys.exit("the kept wisdom was not accepted on reading back")
    from_wisdom = session(lib, FFTW_PATIENT | FFTW_WISDOM_ONLY)
    timings = {n: (estimated[n], from_wisdom[n]) for n in SIZES}
    print(f"{'n':>6s} {'estimated':>12s} {'wisdom':>12s}   (microseconds a transform)")
    for n in SIZES:
        print(f"{n:6d} {estimated[n] * 1e6:12.1f} {from_wisdom[n] * 1e6:12.1f}")

    made_on = f"{platform.system()} {platform.machine()}, {datetime.date.today().isoformat()}"
    with open(args.out, "w") as f:
        f.write(fortran_source(kept[1], version, made_on, timings))
    print(f"wrote {args.out}: {len(kept[1])} characters of wisdom, {version}")


if __name__ == "__main__":
    main()
