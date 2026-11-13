"""
automag.pipeline.readiness
==========================

Polling predicates used by the ``--wait`` mode.

Each stage knows what "the calculations have finished" means for it:

* convergence / linear-response / collinear stages - a result file in
  ``CalcFold`` has accumulated one line per requested calculation;
* relaxation - the ``OUTCAR`` of the relaxation contains VASP's final markers;
* MAE - all ``PhTh_*`` (grid) or ``K_*`` (curve) directories have finished
  ``OUTCAR`` files.

The file-count checks also stop when no progress is observed for
``stable_seconds``.  This avoids waiting forever when a few calculations of a
batch failed, while still allowing the normal post-processing to continue with
the results that are available.

.. codeauthor:: Automag contributors
"""
from __future__ import annotations


import os
import time
from pathlib import Path

from common.results import most_stable_state

FINISH_MARKERS = (
    "General timing and accounting informations",
    "Voluntary context switches",
    "reached required accuracy",
)


def _tail(path, size=16384):
    """Return the last ``size`` bytes of a file (cheap finish-marker check)."""
    with open(path, "rb") as handle:
        handle.seek(0, os.SEEK_END)
        length = handle.tell()
        handle.seek(max(0, length - size))
        return handle.read().decode(errors="ignore")


def outcar_finished(path):
    """Return True if ``path`` looks like a completed VASP OUTCAR."""
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        tail = _tail(path)
    except OSError:
        return False
    return any(marker in tail for marker in FINISH_MARKERS)


def _count_lines(path):
    """Count result lines (lines not starting with whitespace)."""
    path = Path(path)
    if not path.is_file():
        return 0
    try:
        text = path.read_text(errors="ignore")
    except OSError:
        return 0
    return sum(1 for line in text.splitlines()
               if line and not line[0].isspace())


def _wait_for_count(count_fn, expected, *, timeout, interval, stable_seconds,
                    label):
    start = time.monotonic()
    last = -1
    last_change = start
    total = "?" if expected is None else expected
    while True:
        try:
            count = count_fn()
        except OSError:
            count = 0

        if expected is not None and count >= expected:
            print(f"  {label}: {count}/{total} done")
            return count

        now = time.monotonic()
        if count != last:
            last = count
            last_change = now
            print(f"  {label}: {count}/{total} finished")
        elif count > 0 and now - last_change >= stable_seconds:
            print(
                f"WARNING: {label}: {count}/{total} finished and no progress "
                f"for {stable_seconds:.0f}s; continuing with the available results."
            )
            return count

        if now - start >= timeout:
            raise TimeoutError(
                f"Timed out after {timeout:.0f}s waiting for {label} "
                f"({count}/{total})."
            )
        time.sleep(interval)


def _wait_for_predicate(predicate, *, timeout, interval, label):
    start = time.monotonic()
    print(f"  waiting for {label} ...")
    while True:
        try:
            if predicate():
                print(f"  {label}: done")
                return
        except OSError:
            pass
        if time.monotonic() - start >= timeout:
            raise TimeoutError(
                f"Timed out after {timeout:.0f}s waiting for {label}."
            )
        time.sleep(interval)


def wait_for_file(path, *, timeout, interval, label):
    """Wait until ``path`` exists and is non-empty."""
    path = Path(path)
    _wait_for_predicate(
        lambda: path.is_file() and path.stat().st_size > 0,
        timeout=timeout, interval=interval, label=label,
    )


def _conv_expected(config):
    mode = config.get("mode", "kgrid")
    if mode == "encut":
        values = config.get("encut_values")
        return len(list(values)) if values is not None else len(range(500, 1010, 10))
    sigmas = config.get("sigma_values")
    kpts = config.get("kpts_values")
    sigmas = list(sigmas) if sigmas is not None else [0.05, 0.1, 0.15, 0.2]
    kpts = list(kpts) if kpts is not None else list(range(20, 110, 10))
    return len(sigmas) * len(kpts)


def _charges_path(workdir, formula, suffix, calculator):
    base = Path(workdir) / "CalcFold"
    candidates = [
        base / f"{formula}{suffix}_charges_{calculator}.txt",
        base / f"{formula}_charges_{calculator}.txt",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return candidates[0]


def _find_trials(stage_work_dir, workdir, formula, suffix):
    """Locate the ``trials_*`` folder of the running ``2_coll`` instance."""
    names = (f"trials_{formula}{suffix}", f"trials_{formula}")
    roots = []
    if stage_work_dir is not None:
        roots.append(Path(stage_work_dir))
    roots.append(Path(workdir) / "2_coll")
    for root in roots:
        if not root.is_dir():
            continue
        for name in names:
            for trials in root.rglob(name):
                if trials.is_dir():
                    return trials
        for trials in root.rglob("trials_*"):
            if trials.is_dir():
                return trials
    return None


def _count_configurations(stage_work_dir, workdir, formula, suffix):
    trials = _find_trials(stage_work_dir, workdir, formula, suffix)
    if trials is None:
        return None
    total = 0
    for path in trials.glob("configurations*.txt"):
        total += sum(1 for line in path.read_text().splitlines()
                     if line.strip())
    return total


def _mae_configuration(config, workdir, formula, suffix, coll_calculator):
    """Resolve ``configuration`` (possibly ``"lowest"``) for the MAE stage."""
    configuration = str(config["configuration"])
    if configuration.strip().lower() in {"", "lowest", "most_stable",
                                         "most-stable"}:
        return most_stable_state(workdir, formula, suffix, coll_calculator)
    return configuration


def _mae_dirs(workdir, formula, suffix, calculator, configuration):
    """All MAE calculation directories of one configuration.

    A parameter sweep (a list in ``params``) creates one directory per
    parameter combination, so the expected number of calculations has to be
    counted over all of them.
    """
    base = (Path(workdir) / "CalcFold" / f"{formula}{suffix}"
            / calculator / str(configuration))
    if not base.is_dir():
        return []
    return sorted(base.glob("mae_*"))


def _mae_dir(workdir, formula, suffix, calculator, configuration):
    """The most recently modified MAE directory (or ``None``)."""
    candidates = _mae_dirs(workdir, formula, suffix, calculator, configuration)
    if not candidates:
        return None
    return max(candidates, key=lambda item: item.stat().st_mtime)


def _count_finished(parent, pattern):
    if parent is None or not parent.is_dir():
        return 0
    return sum(1 for directory in parent.glob(pattern)
               if outcar_finished(directory / "singlepoint" / "OUTCAR"))


def task_progress(stage, phase, config, workdir, formula, *,
                  stage_work_dir=None):
    """Return ``(finished, expected, label)`` for a phase that submits jobs.

    This is the non-blocking half of :func:`wait_for_task`: the scheduler uses
    it to decide, on every pass, whether the calculations of a step are
    complete.  ``expected`` is ``None`` when the number cannot be derived from
    the settings.
    """
    workdir = Path(workdir)
    calculator = config.get("calculator", "vasp")
    coll_calculator = config.get("coll_calculator") or calculator
    suffix = config.get("struct_suffix", "")

    if stage == "conv_tests" and phase == "submit":
        mode = config.get("mode", "kgrid")
        path = workdir / "CalcFold" / f"{formula}{suffix}_{mode}_{calculator}.txt"
        return _count_lines(path), _conv_expected(config), f"conv_tests {mode}"

    if stage == "lin_response" and phase == "submit":
        path = _charges_path(workdir, formula, suffix, calculator)
        return (_count_lines(path), len(config.get("perturbations", [])),
                "linear response")

    if stage == "coll" and phase == "submit":
        path = (workdir / "CalcFold"
                / f"{formula}{suffix}_singlepoint_{calculator}.txt")
        return (_count_lines(path),
                _count_configurations(stage_work_dir, workdir, formula, suffix),
                "collinear")

    if stage == "relax" and phase == "submit":
        state = config.get("configuration")
        if state and str(state).strip().lower() in {
                "lowest", "most_stable", "most-stable"}:
            state = None
        if not state:
            state = most_stable_state(workdir, formula, suffix, coll_calculator)
        outcar = (workdir / "CalcFold" / f"{formula}{suffix}" / calculator
                  / f"relax_{state}" / "singlepoint" / "OUTCAR")
        return (1 if outcar_finished(outcar) else 0, 1, f"relaxation of {state}")

    if stage == "mae" and phase == "submit":
        configuration = _mae_configuration(config, workdir, formula, suffix,
                                           coll_calculator)
        directories = _mae_dirs(workdir, formula, suffix, calculator,
                                 configuration)
        expected = (int(config["Nph"]) + 1) * (int(config["Nth"]) + 1) \
            * max(1, len(directories))
        finished = sum(_count_finished(directory, "PhTh_*")
                       for directory in directories)
        return finished, expected, "MAE theta/phi"

    if stage == "mae" and phase == "curve":
        configuration = _mae_configuration(config, workdir, formula, suffix,
                                           coll_calculator)
        directories = _mae_dirs(workdir, formula, suffix, calculator,
                                 configuration)
        expected = (int(config["N_MAE"]) + 1) * max(1, len(directories))
        finished = sum(_count_finished(directory, "K_*_RtMAE_*")
                       for directory in directories)
        return finished, expected, "MAE curve"

    return 0, None, f"{stage}.{phase}"


def task_finished(stage, phase, config, workdir, formula, *, stage_work_dir=None):
    """True when the result files show that ``stage.phase`` is complete."""
    finished, expected, _ = task_progress(
        stage, phase, config, workdir, formula, stage_work_dir=stage_work_dir)
    if expected is None:
        return False
    return finished >= expected


def wait_for_task(stage, phase, config, workdir, formula, *, timeout, interval,
                  stable_seconds, stage_work_dir=None):
    """Block until the calculations submitted by ``stage.phase`` are done.

    ``stage_work_dir`` is the directory of the running stage *instance*; it is
    used to find the ``trials_*`` folder when several ``2_coll`` instances
    exist.  Results that live in ``CalcFold`` are always global.
    """
    workdir = Path(workdir)
    calculator = config.get("calculator", "vasp")
    coll_calculator = config.get("coll_calculator") or calculator
    suffix = config.get("struct_suffix", "")

    if stage == "conv_tests" and phase == "submit":
        mode = config.get("mode", "kgrid")
        path = workdir / "CalcFold" / f"{formula}{suffix}_{mode}_{calculator}.txt"
        expected = _conv_expected(config)
        print(f"Waiting for {expected} {mode} convergence calculations ({path}) ...")
        _wait_for_count(lambda: _count_lines(path), expected,
                        timeout=timeout, interval=interval,
                        stable_seconds=stable_seconds,
                        label=f"{stage}.{phase}")

    elif stage == "lin_response" and phase == "submit":
        path = _charges_path(workdir, formula, suffix, calculator)
        expected = len(config.get("perturbations", []))
        print(f"Waiting for {expected} linear-response calculations ({path}) ...")
        _wait_for_count(lambda: _count_lines(path), expected,
                        timeout=timeout, interval=interval,
                        stable_seconds=stable_seconds,
                        label=f"{stage}.{phase}")

    elif stage == "coll" and phase == "submit":
        path = (workdir / "CalcFold"
                / f"{formula}{suffix}_singlepoint_{calculator}.txt")
        expected = _count_configurations(stage_work_dir, workdir, formula, suffix)
        shown = "all" if expected is None else expected
        print(f"Waiting for {shown} collinear calculations ({path}) ...")
        _wait_for_count(lambda: _count_lines(path), expected,
                        timeout=timeout, interval=interval,
                        stable_seconds=stable_seconds,
                        label=f"{stage}.{phase}")

    elif stage == "relax" and phase == "submit":
        state = config.get("configuration")
        if state and str(state).strip().lower() in {
                "lowest", "most_stable", "most-stable"}:
            state = None
        if not state:
            state = most_stable_state(workdir, formula, suffix,
                                      coll_calculator)
        outcar = (workdir / "CalcFold" / f"{formula}{suffix}" / calculator
                  / f"relax_{state}" / "singlepoint" / "OUTCAR")
        _wait_for_predicate(lambda: outcar_finished(outcar),
                            timeout=timeout, interval=interval,
                            label=f"relaxation of {state} ({outcar})")

    elif stage == "mae" and phase == "submit":
        configuration = _mae_configuration(config, workdir, formula, suffix,
                                           coll_calculator)
        directories = _mae_dirs(workdir, formula, suffix, calculator,
                                 configuration)
        expected = (int(config["Nph"]) + 1) * (int(config["Nth"]) + 1) \
            * max(1, len(directories))
        print(f"Waiting for {expected} MAE theta/phi calculations ...")
        _wait_for_count(
            lambda: sum(_count_finished(path, "PhTh_*") for path in directories),
            expected, timeout=timeout, interval=interval,
            stable_seconds=stable_seconds, label=f"{stage}.{phase}")

    elif stage == "mae" and phase == "curve":
        configuration = _mae_configuration(config, workdir, formula, suffix,
                                           coll_calculator)
        directories = _mae_dirs(workdir, formula, suffix, calculator,
                                 configuration)
        expected = (int(config["N_MAE"]) + 1) * max(1, len(directories))
        print(f"Waiting for {expected} MAE curve calculations ...")
        _wait_for_count(
            lambda: sum(_count_finished(path, "K_*_RtMAE_*")
                        for path in directories),
            expected, timeout=timeout, interval=interval,
            stable_seconds=stable_seconds, label=f"{stage}.{phase}")
