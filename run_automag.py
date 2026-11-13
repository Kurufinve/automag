#!/usr/bin/env python
"""
automag.run_automag
===================

Single entry point for the whole pipeline.

Usage::

    python run_automag.py [input_automag.py] [options]

The controlling file may live in any directory and only has to point at the
input structure.  See ``input_automag.py`` for a documented example.

Two internal switches exist so that the same executable can run the stage
scripts and the small helper scripts.  They are used by the frozen
single-file build (``dist/run_automag``) but work in a normal checkout too:

``--run-script <path>``
    execute a stage script (``2_coll/1_submit.py`` and friends).  ``<path>``
    may be relative to the automag root.
``--run-helper <relative path> [args]``
    execute one of the ``common/`` helpers.

.. codeauthor:: Automag contributors
"""

import os
import runpy
import sys


def _bootstrap():
    """Put the automag sources on ``sys.path`` and record where they live."""
    if getattr(sys, "frozen", False):
        root = sys._MEIPASS  # noqa: SLF001 - PyInstaller's unpack directory
    else:
        root = os.path.dirname(os.path.abspath(__file__))
    if root not in sys.path:
        sys.path.insert(0, root)
    # Make the code root explicit: AUTOMAG_PATH may point at a different
    # checkout in the user's shell.
    os.environ["AUTOMAG_CODE_PATH"] = root
    return root


def _take_option(argv, name):
    """Pop ``name <value>`` from ``argv`` and return the value (or None)."""
    if name not in argv:
        return None
    index = argv.index(name)
    try:
        value = argv[index + 1]
    except IndexError:
        raise SystemExit(f"{name} needs a value") from None
    del argv[index:index + 2]
    return value


def main(argv=None):
    argv = list(sys.argv if argv is None else argv)
    root = _bootstrap()

    script = _take_option(argv, "--run-script")
    helper = _take_option(argv, "--run-helper")

    if script or helper:
        relative = script or helper
        path = relative if os.path.isabs(relative) else os.path.join(root, relative)
        if not os.path.isfile(path):
            raise SystemExit(f"Script not found: {path}")
        # The stage scripts read sys.argv themselves (an optional input file).
        sys.argv = [path] + argv[1:]
        runpy.run_path(path, run_name="__main__")
        return 0

    from pipeline.cli import main as cli_main

    return cli_main(argv[1:])


if __name__ == "__main__":
    raise SystemExit(main())
