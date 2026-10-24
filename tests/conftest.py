"""Collection rules for the suites under lucifer/tests.

The keystone suite launches five long jobs from one session fixture. A session fixture
runs once per xdist worker, so a parallel run would launch them several times over into
the same work directories. The suite refuses to start under xdist rather than run twice.
"""

import pytest


def pytest_configure(config: pytest.Config) -> None:
    n = None
    try:
        n = config.getoption("numprocesses")
    except ValueError:
        pass                                   # xdist not installed, nothing to refuse
    if n not in (None, 0, "0"):
        raise pytest.UsageError(
            "xdist is refused for lucifer/tests: the keystone's session fixture runs the "
            "five jobs once, and once per worker would run them into the same work "
            "directories. Run without -n.")
