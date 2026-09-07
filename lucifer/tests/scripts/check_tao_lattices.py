#!/usr/bin/env python3
"""
Every committed lattice loads in Tao.

An FEL segment is a wiggler or undulator carrying Bmad's fel_method attribute, and its
tracking_method stays bmad_standard, so a Lucifer lattice is an ordinary Bmad lattice that
any Bmad program reads. This check holds that: it runs Tao over every committed lattice
and requires a zero exit with no error line. A lattice Tao cannot load fails the benchmark.

The tree under tests/bmad and examples is copied into the work directory with its layout
intact, since a wrapper reaches its base by relative path, and so that no digested file
lands in the source tree.

Usage:
  check_tao_lattices.py --tao <tao> --latdir <lucifer/tests/bmad> --examples <lucifer/examples>
                        --workdir <dir>
"""

import argparse
import pathlib
import re
import shutil
import subprocess
import sys

# A wrapper sets attributes on a line it calls in and names no line of its own, so Tao is
# pointed at the lattices that name one.
RE_USE = re.compile(r'^\s*use\s*,', re.MULTILINE | re.IGNORECASE)

# Written by a run rather than by hand, so outside the committed set.
SKIP_DIRS = {'output'}

# A fragment by design: its own header says the wrappers define ANG before calling it, so
# it names a line but is not a lattice on its own. Its two wrappers are both checked.
SKIP_FILES = {'chicane_base.bmad'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tao', required=True)
    ap.add_argument('--latdir', required=True)
    ap.add_argument('--examples', required=True)
    ap.add_argument('--workdir', required=True)
    args = ap.parse_args()

    wd = pathlib.Path(args.workdir)
    if wd.exists():
        shutil.rmtree(wd)
    wd.mkdir(parents=True)

    # Copy with the layout intact so a relative call resolves as it does in the tree. The
    # wake tables and the reference lattices come too, since an element names one by path.
    roots = {'bmad': pathlib.Path(args.latdir), 'examples': pathlib.Path(args.examples)}
    lats = []
    for tag, root in roots.items():
        for p in sorted(root.rglob('*')):
            if not p.is_file() or SKIP_DIRS & set(p.parts):
                continue
            if p.suffix not in ('.bmad', '.wake', '.lat'):
                continue
            dest = wd / tag / p.relative_to(root)
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(p, dest)
            if p.suffix == '.bmad' and p.name not in SKIP_FILES and RE_USE.search(p.read_text()):
                lats.append(dest)

    if not lats:
        print('FAIL: no lattices found to load', file=sys.stderr)
        return 1

    n_bad = 0
    for p in lats:
        r = subprocess.run([args.tao, '-lat', p.name, '-noplot', '-command', 'exit'],
                           cwd=p.parent, stdin=subprocess.DEVNULL,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        bad = [ln for ln in r.stdout.splitlines() if 'ERROR' in ln.upper()]
        if r.returncode != 0 or bad:
            n_bad += 1
            print(f'  FAIL: {p.name}: Tao exited {r.returncode}')
            for ln in (bad or r.stdout.splitlines())[:4]:
                print(f'         {ln.strip()[:110]}')

    print(f'  Tao loaded {len(lats) - n_bad} of {len(lats)} committed lattices')
    print(f'checks: {"PASS" if n_bad == 0 else "FAIL"}')
    return 1 if n_bad else 0


if __name__ == '__main__':
    sys.exit(main())
