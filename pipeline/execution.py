"""
automag.pipeline.execution
==========================

How the work of the pipeline is started, supervised and resumed.

Two kinds of work exist:

**step scripts**
    every phase (``submit``, ``postprocess``, ``axes``, ...) is implemented by
    one or more automag scripts.  They can be executed

    ``inprocess``   inside the automag process (``runpy``);
    ``process``     as a separate process on the step's host (local or ``ssh``);
    ``scheduler``   as a batch job of the host's task manager (SLURM/PBS), so
                    a long analysis does not occupy the login node.

    ``auto`` (the default) picks ``scheduler`` when the step's ``header``
    contains SLURM/PBS directives, otherwise ``process`` or ``inprocess``
    depending on the resources of the machine.

**calculation jobs**
    a ``submit`` phase writes the inputs of many independent calculations and
    one job script per calculation.  Automag reads that list from the
    *manifest* (:data:`MANIFEST_NAME`, written by :class:`common.SubmitManual`),
    starts at most ``max_parallel`` of them at a time and remembers every
    started job in :class:`pipeline.state.StateStore`.  That is what makes

    * parallel execution of a step (``max_parallel``),
    * priorities (``depends_on``: the 4_mae reference calculation first),
    * and a resumable run (only the calculations without a result are
      restarted) possible.

.. codeauthor:: Automag contributors
"""
from __future__ import annotations


import json
import os
import runpy
import shlex
import subprocess
import sys
import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path

from pipeline.hosts import Host, HostError

#: Name of the manifest written by ``SubmitManual`` in a stage directory.
MANIFEST_NAME = "automag_tasks.json"

#: Environment variable through which the runner tells ``SubmitManual`` to
#: record the job scripts it creates instead of launching them.
MANIFEST_ENV = "AUTOMAG_TASK_MANIFEST"

#: Environment variable carrying the step key into the stage script.
STEP_KEY_ENV = "AUTOMAG_STEP_KEY"

#: Execution modes.
MODES = ("auto", "inprocess", "process", "scheduler")

LOG_DIR = ".automag_logs"
BATCH_DIR = ".automag_batches"

#: Marker written by a wrapper batch script once the step script succeeded.
STEP_DONE_SUFFIX = ".done"


# ---------------------------------------------------------------------------
# Manifest (SubmitManual -> runner)
# ---------------------------------------------------------------------------
def manifest_path(work_dir):
    return Path(work_dir) / MANIFEST_NAME


def read_manifest(work_dir):
    """Return the manifest entries of a stage directory (possibly empty)."""
    path = manifest_path(work_dir)
    if not path.is_file():
        return []
    entries = []
    for number, line in enumerate(path.read_text().splitlines(), start=1):
        line = line.strip()
        if not line:
            continue
        try:
            entries.append(json.loads(line))
        except ValueError:
            print(f"WARNING: {path}:{number}: ignoring a malformed line")
    return entries


def write_manifest(work_dir, entries):
    """(Re)write the manifest of a stage directory."""
    path = manifest_path(work_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as handle:
        for entry in entries:
            handle.write(json.dumps(entry, default=str) + "\n")
    return path


def clear_manifest(work_dir):
    path = manifest_path(work_dir)
    if path.is_file():
        path.unlink()


# ---------------------------------------------------------------------------
# Mode selection
# ---------------------------------------------------------------------------
def scheduler_directives(header):
    """Return ``"slurm"`` / ``"pbs"`` when the batch header asks for it."""
    if not header:
        return None
    text = str(header)
    lowered = text.lower()
    if "#sbatch" in lowered or "sbatch " in lowered:
        return "slurm"
    if "#pbs" in lowered or "qsub " in lowered:
        return "pbs"
    return None


def machine_resources():
    """CPU count and available memory of this machine (best effort)."""
    cpus = os.cpu_count() or 1
    try:
        cpus = len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        pass
    available = None
    try:
        with open("/proc/meminfo") as handle:
            for line in handle:
                if line.startswith("MemAvailable:"):
                    available = int(line.split()[1]) * 1024
                    break
    except OSError:
        available = None
    return {"cpus": cpus, "memory_bytes": available}


def can_spawn_process(min_cpus=2, min_memory=512 * 1024 * 1024):
    """True when a separate process is affordable on this machine."""
    resources = machine_resources()
    if resources["cpus"] < min_cpus:
        return False
    memory = resources["memory_bytes"]
    return memory is None or memory >= min_memory


def resolve_mode(instance, phase, host=None, override=None):
    """Return how the *step script itself* is executed.

    ``auto`` selects ``process`` (a separate automag process) or - on a machine
    that cannot afford another process - ``inprocess`` (inside the automag
    process).  ``scheduler`` (the script becomes a batch job) is never chosen
    automatically: what the batch directives in the step ``header`` select is
    how the *calculations* of the step are started (see :func:`launcher_for`).
    """
    plan = instance.phase_execution(phase) if instance is not None else {}
    mode = str(override or plan.get("mode") or "auto").lower()
    if mode not in MODES:
        raise HostError(
            f"Unknown execution mode '{mode}' for phase '{phase}'. "
            f"Use one of {', '.join(MODES)}.")
    if mode != "auto":
        return mode
    if host is not None and not host.is_local:
        return "process"
    return "process" if can_spawn_process() else "inprocess"


def _header_of(instance, phase=None):
    if instance is None:
        return ""
    plan = instance.phase_execution(phase) if phase else {}
    return plan.get("header") or instance.jobheader


def launcher_for(entry, host, header=""):
    """How one calculation job is started: ``sbatch``/``qsub``/``process``."""
    launcher = str(entry.get("launcher") or "auto").lower()
    if launcher in ("sbatch", "qsub", "process", "sh"):
        return launcher
    if host is not None and host.scheduler in ("slurm", "pbs"):
        return "sbatch" if host.scheduler == "slurm" else "qsub"
    directives = scheduler_directives(header)
    if directives == "slurm":
        return "sbatch"
    if directives == "pbs":
        return "qsub"
    return "process"


# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------
def log_path(work_dir, step_key):
    directory = Path(work_dir) / LOG_DIR
    directory.mkdir(parents=True, exist_ok=True)
    return directory / f"{step_key.replace('.', '_')}.log"


def step_sentinel(work_dir, step_key):
    """File a wrapper batch script touches when the step finished."""
    return log_path(work_dir, step_key).with_suffix(".done")


def append_log(path, text):
    if path is None or not text:
        return
    try:
        with open(path, "a") as handle:
            handle.write(text if text.endswith("\n") else text + "\n")
    except OSError:
        pass


# ---------------------------------------------------------------------------
# Running step scripts
# ---------------------------------------------------------------------------
def run_inprocess(script_path, work_dir, env=None, argv=None):
    """Execute ``script_path`` inside this process, like ``python script.py``.

    Returns the exit code (0 on success).  Exceptions of the script are logged
    and reported as a failure instead of killing automag.
    """
    script_path = Path(script_path)
    work_dir = Path(work_dir)
    previous_cwd = os.getcwd()
    previous_env = {}
    for key, value in (env or {}).items():
        if key in os.environ:
            previous_env[key] = os.environ[key]
        else:
            previous_env[key] = None
        os.environ[key] = str(value)
    sys.path.insert(0, str(work_dir))
    previous_argv = sys.argv
    sys.argv = [str(script_path)] + list(argv or [])
    try:
        os.chdir(str(work_dir))
        runpy.run_path(str(script_path), run_name="__main__")
        return 0
    except SystemExit as exit_signal:
        code = exit_signal.code
        return 0 if code in (None, 0) else int(code) if isinstance(code, int) else 1
    except BaseException:  # noqa: BLE001 - a broken step must not kill automag
        traceback.print_exc()
        return 1
    finally:
        sys.argv = previous_argv
        try:
            os.chdir(previous_cwd)
        except OSError:
            pass
        if str(work_dir) in sys.path:
            sys.path.remove(str(work_dir))
        for key, value in previous_env.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def run_sync(host, command, cwd, log=None, timeout=None, env=None):
    """Run ``command`` on ``host`` in the foreground; return ``(rc, out, err)``."""
    line = host.wrap(command, cwd=cwd, env=env)
    try:
        completed = subprocess.run(
            line, shell=True, text=True, capture_output=True, timeout=timeout,
            executable="/bin/bash" if os.path.exists("/bin/bash") else None)
    except subprocess.TimeoutExpired as error:
        append_log(log, f"TIMEOUT after {timeout}s: {line}")
        return 124, "", str(error)
    if log is not None:
        append_log(log, f"\n=== {time.strftime('%Y-%m-%d %H:%M:%S')} @{host.name}\n"
                        f"$ {line}\n{completed.stdout or ''}{completed.stderr or ''}")
        if completed.returncode != 0:
            append_log(log, f"exit code: {completed.returncode}")
    if completed.stdout:
        sys.stdout.write(completed.stdout)
        if not completed.stdout.endswith("\n"):
            sys.stdout.write("\n")
    if completed.returncode != 0 and completed.stderr:
        sys.stderr.write(completed.stderr)
    return completed.returncode, completed.stdout or "", completed.stderr or ""


def automag_command(host, relative_script):
    """Command that runs an automag stage script on ``host``."""
    if getattr(sys, "frozen", False) and host.is_local:
        return f"{shlex.quote(sys.executable)} --run-script {relative_script}"
    return host.automag_command(relative_script)


# ---------------------------------------------------------------------------
# Completion probes
# ---------------------------------------------------------------------------
def entry_markers(entry):
    """Return the list of ``done`` specifications of a manifest entry."""
    done = entry.get("done") if isinstance(entry, dict) else entry
    if isinstance(done, str):
        return [{"kind": "file", "path": done}]
    if isinstance(done, dict):
        return [done] if done.get("path") else []
    if isinstance(done, (list, tuple)):
        return [item for item in done
                if isinstance(item, dict) and item.get("path")]
    return []


def marker_mtime(path):
    try:
        return Path(path).stat().st_mtime
    except OSError:
        return None


def _spec_satisfied(done, min_mtime=None, base=None):
    from pipeline import readiness

    path = Path(done["path"])
    if base is not None and not path.is_absolute():
        # 'singlepoint/OUTCAR' style marker relative to the calculation dir
        path = Path(base) / path
    kind = done.get("kind", "file")
    if kind == "outcar":
        satisfied = readiness.outcar_finished(path)
    else:
        # existence is the signal: the job writes the sentinel with `touch`
        # (zero bytes) and the result file as its very last action
        try:
            satisfied = path.is_file()
        except OSError:
            satisfied = False
    if not satisfied:
        return False
    if min_mtime is None:
        return True
    stamp = marker_mtime(path)
    return stamp is not None and stamp >= min_mtime - 1e-3


def marker_satisfied(entry, min_mtime=None):
    """True when the calculation of ``entry`` produced its result marker.

    ``min_mtime`` ignores markers older than a restart of the step, so that
    changing a parameter really recomputes the stale calculation.
    """
    specs = entry_markers(entry)
    if not specs:
        return False
    base = entry.get("jobdir")
    return any(_spec_satisfied(spec, min_mtime, base=base) for spec in specs)


# ---------------------------------------------------------------------------
# Calculation jobs of one step
# ---------------------------------------------------------------------------
@dataclass
class Job:
    """One calculation of a ``submit`` phase."""

    jobdir: str
    jobscript: str
    launcher: str = "auto"
    kind: str = "calculation"
    depends_on: str | None = None
    done: dict | None = None
    mode: str = ""
    raw: dict = field(default_factory=dict)

    @property
    def entry(self):
        return dict(self.raw or {
            "jobdir": self.jobdir, "jobscript": self.jobscript,
            "launcher": self.launcher, "kind": self.kind,
            "depends_on": self.depends_on, "done": self.done,
        })

    @classmethod
    def from_entry(cls, entry):
        return cls(
            jobdir=str(entry.get("jobdir")),
            jobscript=str(entry.get("jobscript") or "jobscript"),
            launcher=str(entry.get("launcher") or "auto"),
            kind=str(entry.get("kind") or "calculation"),
            depends_on=(str(entry["depends_on"])
                        if entry.get("depends_on") else None),
            done=entry.get("done"),
            mode=str(entry.get("mode") or ""),
            raw=dict(entry),
        )


def manifest_jobs(work_dir):
    """Read the manifest as :class:`Job` objects, dropping duplicates."""
    jobs, seen = [], set()
    for entry in read_manifest(work_dir):
        if not entry.get("jobdir"):
            continue
        key = (str(entry["jobdir"]), str(entry.get("jobscript")))
        if key in seen:
            continue
        seen.add(key)
        jobs.append(Job.from_entry(entry))
    return jobs


class StepJobs:
    """Submit, supervise and resume the calculations of one ``submit`` phase.

    Every calculation is a separate job, so ``max_parallel`` really is the
    number of concurrent calculations.  A job whose ``depends_on`` is not
    finished yet is never started (that is how the 4_mae reference calculation
    is ordered before the angle grid).
    """

    def __init__(self, state, step_key, host, work_dir, header="",
                 calculator=None, max_parallel=1, max_attempts=3,
                 dry_run=False, min_mtime=None, log=None):
        self.state = state
        self.step_key = step_key
        self.host = host
        self.work_dir = Path(work_dir)
        self.header = header or ""
        self.calculator = calculator
        self.max_parallel = max(1, int(max_parallel or 1))
        self.max_attempts = max(1, int(max_attempts or 3))
        self.dry_run = dry_run
        self.min_mtime = min_mtime
        self.log = log
        self.jobs = {}

    # -- registration ---------------------------------------------------
    def register(self, jobs=None):
        """Register the manifest calculations as items of the step."""
        if jobs is None:
            jobs = manifest_jobs(self.work_dir)
        self.jobs = {job.jobdir: job for job in jobs}
        paths = [job.jobdir for job in jobs]
        added = self.state.add_items(self.step_key, paths)
        if added:
            print(f"    {len(added)} calculation(s) registered")
        return added

    def count(self):
        return len(self.jobs)

    # -- bookkeeping ----------------------------------------------------
    def _task_items(self):
        items = set()
        for task in self.state.tasks(self.step_key,
                                     status=("pending", "submitted", "running")):
            items.update(task.get("items") or [])
        return items

    def dependencies_finished(self, job, done_paths):
        if not job.depends_on:
            return True
        if job.depends_on in done_paths:
            return True
        target = self.jobs.get(job.depends_on)
        if target is not None:
            return marker_satisfied(target.entry, self.min_mtime)
        dependency_marker = (job.raw or {}).get("dep_done")
        if dependency_marker:
            return marker_satisfied({"done": dependency_marker}, self.min_mtime)
        return False

    def refresh(self):
        """Mark finished calculations as done; return the pending items."""
        pending = self.state.pending_items(self.step_key)
        still = []
        for item in pending:
            job = self.jobs.get(item["path"])
            if job is not None and marker_satisfied(job.entry, self.min_mtime):
                self.state.set_item(item["item_key"], status="done")
                continue
            still.append(item)
        return still

    def progress(self):
        return self.state.items_progress(self.step_key)

    def running_tasks(self):
        return self.state.tasks(self.step_key,
                               status=("pending", "submitted", "running"))

    def sync_tasks(self):
        """Close the tasks whose calculations finished or whose job died.

        The item markers are the source of truth: a task is finished as soon
        as every calculation it started has a result.  Jobs the host no longer
        knows (process gone, job purged from the queue) are closed as well, so
        that a missing result is retried instead of blocking a slot forever.
        """
        for task in self.running_tasks():
            items = [self.state.get_item(key) for key in task.get("items") or []]
            if items and all(item is not None and item["status"] == "done"
                             for item in items):
                self.state.update_task(task["task_key"], status="finished")
                continue
            info = self.host.job_state(task["external_id"], kind=task["kind"],
                                       dry_run=self.dry_run)
            state = info.get("state")
            if state == "failed":
                self.state.update_task(task["task_key"], status="failed",
                                       detail={"output": info.get("detail")})
                print(f"    WARNING: job {task['external_id']} failed")
            elif state in ("finished", "unknown"):
                self.state.update_task(task["task_key"], status="finished",
                                       detail={"output": info.get("detail")})
                for key in task.get("items") or []:
                    item = self.state.get_item(key)
                    if item is not None and item["status"] == "running":
                        self.state.set_item(key, status="pending")

    # -- submission -----------------------------------------------------
    def pump(self, pending=None):
        """Start as many pending calculations as ``max_parallel`` allows."""
        self.sync_tasks()
        if pending is None:
            pending = self.refresh()
        if not pending:
            return 0
        busy_items = self._task_items()
        free = self.max_parallel - len(self.running_tasks())
        if free <= 0:
            return 0
        done_paths = {item["path"] for item in self.state.items(self.step_key)
                      if item["status"] == "done"}
        started = 0
        for item in pending:
            if started >= free:
                break
            if item["item_key"] in busy_items:
                continue
            if item["attempts"] >= self.max_attempts:
                continue
            job = self.jobs.get(item["path"])
            if job is None:
                print(f"    WARNING: {item['path']} is not in the manifest; "
                      "skipping it")
                self.state.set_item(item["item_key"], status="failed")
                continue
            if not self.dependencies_finished(job, done_paths):
                continue
            if self.submit_job(job, item):
                started += 1
        return started

    def submit_job(self, job, item):
        """Start one calculation and record it in the state."""
        launcher = launcher_for(job.entry, self.host, self.header)
        script = str(Path(job.jobdir) / job.jobscript)
        kind = ("slurm" if launcher == "sbatch"
                else "pbs" if launcher == "qsub" else "process")
        command = (self.host.submit_command(script, cwd=job.jobdir,
                                            launcher=launcher)
                   if kind != "process"
                   else f"cd {shlex.quote(job.jobdir)} && sh {shlex.quote(job.jobscript)}")
        task_key = self.state.add_task(
            self.step_key, host=self.host.name, kind=kind, workdir=job.jobdir,
            command=command, items=[item["item_key"]], launcher=launcher)
        self.state.bump_item_attempts([item["item_key"]], task_key=task_key)
        label = Path(job.jobdir).name
        if self.dry_run:
            print(f"    [dry-run @{self.host.name}] {kind}: {command}")
            self.state.update_task(task_key, status="finished")
            self.state.set_item(item["item_key"], status="done")
            return True
        try:
            if kind == "process":
                pid = self.host.start_background(
                    f"sh {shlex.quote(job.jobscript)}", cwd=job.jobdir,
                    log_path=str(log_path(self.work_dir, self.step_key)))
                external_id = pid
                if not pid:
                    external_id = None
            else:
                external_id, output = self.host.submit(
                    script, cwd=job.jobdir, launcher=launcher)
        except (HostError, OSError, subprocess.SubprocessError) as error:
            print(f"    WARNING: could not start {job.jobdir}: {error}")
            self.state.update_task(task_key, status="failed",
                                   detail={"error": str(error)})
            return False
        self.state.update_task(
            task_key, status="running", external_id=external_id,
            detail={"started_at": time.time(),
                    "log": str(log_path(self.work_dir, self.step_key))})
        print(f"    started {label} ({kind}"
              + (f", id {external_id}" if external_id else "") + ")")
        return True

    # -- supervision ----------------------------------------------------
    def poll(self):
        """One polling round; return True when nothing is left to compute."""
        self.sync_tasks()
        pending = self.refresh()
        done, total = self.progress()
        if total:
            print(f"    {done}/{total} calculation(s) done")
        if not pending:
            return True
        if self.state.tasks(self.step_key,
                            status=("pending", "submitted", "running")):
            self.pump(pending)
            return False
        retry = [item for item in pending
                 if item["attempts"] < self.max_attempts]
        exhausted = [item for item in pending
                     if item["attempts"] >= self.max_attempts]
        if exhausted:
            print(f"    WARNING: {len(exhausted)} calculation(s) still missing "
                  f"after {self.max_attempts} attempt(s):")
            for item in exhausted[:5]:
                print(f"      - {Path(item['path']).name}")
            if len(exhausted) > 5:
                print(f"      ... and {len(exhausted) - 5} more")
        if retry:
            print(f"    resubmitting {len(retry)} missing calculation(s)")
            self.pump(retry)
            return False
        return True

    def wait(self, timeout=None, interval=30, stable_seconds=None,
             expected=None):
        """Poll until every calculation is done or out of attempts."""
        started = time.monotonic()
        last_progress, last_done = started, -1
        while True:
            if self.poll():
                return True
            done, total = self.progress()
            moment = time.monotonic()
            if done != last_done:
                last_done, last_progress = done, moment
            if (stable_seconds is not None and done > 0
                    and done == last_done
                    and moment - last_progress >= stable_seconds
                    and not self.state.tasks(
                        self.step_key, status=("pending", "submitted",
                                               "running"))):
                print(f"    WARNING: no progress for {stable_seconds:.0f}s; "
                      "continuing with the available results")
                return True
            if timeout is not None and moment - started >= timeout:
                raise TimeoutError(
                    f"Timed out waiting for {self.step_key} "
                    f"({done}/{total} done).")
            time.sleep(interval)


# ---------------------------------------------------------------------------
# Step scripts as batch jobs
# ---------------------------------------------------------------------------
def step_wrapper_script(host, work_dir, step_key, command, header=None):
    """Write a batch script that runs one automag step script."""
    directory = Path(work_dir) / BATCH_DIR
    directory.mkdir(parents=True, exist_ok=True)
    script = directory / f"{step_key.replace('.', '_')}.sh"
    sentinel = step_sentinel(work_dir, step_key)
    lines = []
    if header:
        lines.append(str(header).rstrip("\n"))
    else:
        lines.append("#!/bin/sh")
    lines.append(f"# automag step {step_key}")
    lines.append(f"rm -f {shlex.quote(str(sentinel))}")
    lines.append("cd " + shlex.quote(str(work_dir)))
    lines.extend(host.shell_prefix())
    lines.append(command)
    lines.append(f"status=$?")
    lines.append(f'[ "$status" -eq 0 ] && touch {shlex.quote(str(sentinel))}')
    lines.append('exit "$status"')
    script.write_text("\n".join(lines) + "\n")
    try:
        os.chmod(script, 0o755)
    except OSError:
        pass
    return script, sentinel
