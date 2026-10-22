# Lucifer

Lucifer is a free-electron-laser tracker built on Bmad. Its physics is transcribed from Genesis 1.3 Version 4 and validated against it from bitwise-identical starting states. The program, its documentation and its checks are in [`lucifer/`](lucifer/README.md), and the documentation is published at <https://christophermayes.github.io/lucifer/>.

Lucifer is under development and unreleased. It has no version number, and its inputs, outputs and interfaces may change without notice.

Lucifer is a Bmad user program and builds against a Bmad checkout at or after the release named in [`lucifer/BMAD_RELEASE`](lucifer/BMAD_RELEASE), which `BMAD_DIR` names. The build and the checks are described in [`lucifer/doc/user-guide.md`](lucifer/doc/user-guide.md#building) and [`lucifer/doc/validation.md`](lucifer/doc/validation.md).

## License

Copyright 2026 the Lucifer authors.

Lucifer is licensed under the GNU General Public License version 3 only (SPDX: `GPL-3.0-only`). The license text is in [`LICENSE`](LICENSE).

The routines transcribed from Genesis 1.3 Version 4 are derivative works of it, modified from 2026-07-31 on, and are distributed under the same license. Genesis 1.3 Version 4 carries the notice

```
Copyright 2010, 2011 Paul Scherrer Institut. All rights reserved.
```

and is licensed under the GNU General Public License version 3. Every transcribed routine cites the Genesis4 file and lines it came from, and those citations are kept.
