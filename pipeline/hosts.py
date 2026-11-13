"""
automag.pipeline.hosts
======================

Where a pipeline step runs.

A *host* describes one machine that can execute automag commands and submit
jobs.  Each stage instance picks a host by name (``"host": "cluster1"``) or
defines one inline.  Hosts may be mixed freely in one pipeline: one step
locally, the next on cluster A, the next on cluster B.

Configuration
-------------
Named hosts live in the pipeline input file::

    PIPELINE = {
        "hosts": {
            "local":   {"transport": "local"},
            "cluster": {"transport": "ssh",
                        "host": "login.cluster.example",
                        "user": "me",
                        "identity_file": "~/.ssh/id_ed25519",
                        "automag_home": "/shared/automag",
                        "python": "python3",
                        "scheduler": "slurm"},          # slurm | pbs | none
            "gpu":     {"transport": "shell",           # commands through a wrapper
                        "wrapper": "docker exec -i gpuc"},
        },
        "run": [
            {"id": "coll", "stage": "coll", "host": "cluster", "params": {...}},
            {"id": "mae",  "stage": "mae",
             "host": {"transport": "ssh", "host": "other.example", "user": "me",
                      "automag_home": "/opt/automag"},
             "params": {...}},
            {"id": "relax", "stage": "relax", "host": "local"},
        ],
    }

Transports
----------
``local``
    subprocess on this machine (the default).
``ssh``
    the same commands prefixed with ``ssh``; requires the *project directory*
    (``workdir``) and an automag checkout/executable to be reachable from the
    remote side.
``shell``
    commands are run through ``wrapper`` on this machine, e.g. inside a
    container.  Mainly useful for testing an ssh-like setup without a second
    host.

``run_automag.py --list-hosts`` prints the resolved table.

.. codeauthor:: Automag contributors
"""
from __future__ import annotations


import os
import shlex
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

TRANSPORTS = ("local", "ssh", "shell")
SCHEDULERS = ("slurm", "pbs", "none")

#: Default location of the automag sources / executable on a remote host.
DEFAULT_AUTOMAG_HOME = "~/automag"

SSH_OPTIONS = ("-o", "BatchMode=yes")


class HostError(RuntimeError):
    pass


@dataclass
class Host:
    """One machine that can run automag commands and submit jobs."""

    name: str = "local"
    transport: str = "local"
    host: str | None = None
    user: str | None = None
    port: int | None = None
    identity_file: str | None = None
    workdir: str | None = None
    """Project directory on that machine (defaults to the local workdir)."""
    automag_home: str | None = None
    automag_exe: str | None = None
    python: str = "python3"
    scheduler: str | None = None
    environment: tuple = ()
    """Shell lines executed before a command (``module load``, ``source``...)."""
    wrapper: str | None = None
    ssh_options: tuple = ()
    extra: dict = field(default_factory=dict)

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------
    @classmethod
    def from_config(cls, name, config):
        if config is None:
            config = {}
        if isinstance(config, str):
            config = {"transport": config}
        if not isinstance(config, dict):
            raise HostError(f"Host '{name}' must be a mapping or a string.")
        values = dict(config)
        transport = str(values.pop("transport", "local")).lower()
        if transport not in TRANSPORTS:
            raise HostError(
                f"Host '{name}': unknown transport '{transport}'. "
                f"Use one of {', '.join(TRANSPORTS)}.")
        options = values.pop("ssh_options", ()) or ()
        if isinstance(options, str):
            options = tuple(shlex.split(options))
        environment = values.pop("environment", ()) or ()
        if isinstance(environment, str):
            environment = (environment,)
        scheduler = values.pop("scheduler", None)
        if scheduler is not None:
            scheduler = str(scheduler).lower()
            if scheduler not in SCHEDULERS:
                raise HostError(
                    f"Host '{name}': unknown scheduler '{scheduler}'. "
                    f"Use one of {', '.join(SCHEDULERS)}.")
        port = values.pop("port", None)
        known = {
            "host", "user", "identity_file", "workdir", "automag_home",
            "automag_exe", "python", "wrapper",
        }
        kwargs = {key: values.pop(key) for key in list(values) if key in known}
        return cls(
            name=name, transport=transport, scheduler=scheduler,
            port=int(port) if port else None,
            ssh_options=tuple(options), environment=tuple(environment),
            extra=values, **kwargs,
        )

    def resolved(self, local_workdir=None):
        """Fill in the paths that default to the local ones."""
        workdir = self.workdir or (str(local_workdir) if local_workdir else None)
        return Host(**{**self.__dict__, "workdir": workdir,
                       "environment": tuple(self.environment),
                       "ssh_options": tuple(self.ssh_options),
                       "extra": dict(self.extra)}, )

    # ------------------------------------------------------------------
    # Introspection
    # ------------------------------------------------------------------
    @property
    def is_local(self):
        return self.transport == "local"

    def describe(self):
        if self.transport == "local":
            target = "this machine"
        elif self.transport == "ssh":
            target = f"{self.user + '@' if self.user else ''}{self.host}"
            if self.port:
                target += f":{self.port}"
        else:
            target = f"wrapper: {self.wrapper}"
        scheduler = self.scheduler or "none"
        return f"{self.name}: {self.transport} -> {target} (scheduler: {scheduler})"

    # ------------------------------------------------------------------
    # Remote automag invocation
    # ------------------------------------------------------------------
    def automag_command(self, relative_script, extra_args=""):
        """Command that runs an automag stage script on this host."""
        if self.automag_exe:
            command = f"{self.automag_exe} --run-script {relative_script}"
        elif self.transport in ("local", "shell") and (
                self.automag_home is None or
                self.python in (None, "", "python3", "python")):
            # Use the interpreter and the checkout that are running automag
            # right now (a bare 'python3' may be a different environment).
            import sys

            from common.paths import code_root

            home = self.automag_home or str(code_root())
            interpreter = (self.python if self.python not in
                           (None, "", "python3", "python") else sys.executable)
            command = (f"{shlex.quote(interpreter)} "
                       f"{shlex.quote(home.rstrip('/') + '/run_automag.py')} "
                       f"--run-script {relative_script}")
        else:
            home = self.automag_home or DEFAULT_AUTOMAG_HOME
            command = (f"{self.python} {home.rstrip('/')}/run_automag.py "
                       f"--run-script {relative_script}")
        if extra_args:
            command = f"{command} {extra_args}"
        return command

    def shell_prefix(self):
        """Environment lines that must precede every command."""
        return [line for line in self.environment if line.strip()]

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------
    def wrap(self, command, cwd=None, env=None):
        """Return the full command line to run ``command`` on this host."""
        lines = []
        directory = cwd if cwd is not None else self.workdir
        if directory:
            lines.append(f"cd {shlex.quote(str(directory))}")
        lines.extend(self.shell_prefix())
        if env:
            for key, value in env.items():
                if value is None:
                    continue
                lines.append(f"export {key}={shlex.quote(str(value))}")
        lines.append(command)
        joined = " && ".join(lines)

        if self.transport == "local":
            return joined
        if self.transport == "shell":
            if not self.wrapper:
                raise HostError(
                    f"Host '{self.name}': transport 'shell' needs a 'wrapper' "
                    "(for example 'docker exec -i container').")
            return f"{self.wrapper} sh -c {shlex.quote(joined)}"
        # ssh
        if not self.host:
            raise HostError(f"Host '{self.name}': transport 'ssh' needs 'host'.")
        target = f"{self.user}@{self.host}" if self.user else self.host
        options = list(SSH_OPTIONS) + list(self.ssh_options)
        if self.port:
            options += ["-p", str(self.port)]
        if self.identity_file:
            options += ["-i", os.path.expanduser(str(self.identity_file))]
        option_string = " ".join(shlex.quote(item) for item in options)
        return f"ssh {option_string} {shlex.quote(target)} sh -c {shlex.quote(joined)}"

    def run(self, command, cwd=None, timeout=None, check=False, capture=True,
            dry_run=False):
        """Run ``command`` on this host."""
        line = self.wrap(command, cwd=cwd)
        if dry_run:
            print(f"    [dry-run @{self.name}] {line}")
            return _completed(0, "", "")
        if self.transport == "local":
            return subprocess.run(
                line, shell=True, timeout=timeout, check=check,
                capture_output=capture, text=True,
                executable="/bin/bash" if os.path.exists("/bin/bash") else None)
        return subprocess.run(
            line, shell=True, timeout=timeout, check=check,
            capture_output=capture, text=True)

    def submit_command(self, script, cwd=None, launcher=None):
        """Scheduler submission line for a job script on this host."""
        scheduler = self.scheduler or "none"
        if launcher:
            launcher = str(launcher)
        elif scheduler == "slurm":
            launcher = "sbatch"
        elif scheduler == "pbs":
            launcher = "qsub"
        else:
            launcher = "sh"
        directory = cwd if cwd is not None else self.workdir
        prefix = f"cd {shlex.quote(str(directory))} && " if directory else ""
        return f"{prefix}{launcher} {shlex.quote(str(script))}"

    def submit(self, script, cwd=None, timeout=120, dry_run=False, launcher=None):
        """Submit a job script and return ``(external_id, output)``."""
        line = self.submit_command(script, cwd=cwd, launcher=launcher)
        result = self.run(line, timeout=timeout, dry_run=dry_run)
        text = (result.stdout or "") + (result.stderr or "")
        if result.returncode != 0:
            raise HostError(
                f"Host '{self.name}': submission failed for {script}:\n{text}")
        job_id = parse_job_id(text, self.scheduler)
        if job_id is None and launcher in (None, "sh", "bash"):
            # A plain 'sh jobscript' runs synchronously: the job is already
            # finished, so report the launching pid-less completion.
            return None, text
        return job_id, text

    # ------------------------------------------------------------------
    # Supervision of already started work
    # ------------------------------------------------------------------
    def job_state(self, external_id, kind=None, timeout=30, dry_run=False):
        """Return the state of a previously submitted job.

        ``{"state": "running" | "pending" | "finished" | "failed" | "unknown"}``
        ``"unknown"`` means the scheduler no longer knows the job; the caller
        then falls back to the calculation's own result marker.
        """
        if dry_run or external_id in (None, ""):
            return {"state": "unknown", "detail": "not submitted"}
        scheduler = (kind or self.scheduler or "none").lower()
        if scheduler in ("none", "process", "pid", ""):
            command = f"ps -p {shlex.quote(str(external_id))} -o stat="
        elif scheduler == "slurm":
            command = f"squeue -h -j {shlex.quote(str(external_id))} -o %T"
        else:  # pbs / pbspro / torque
            command = f"qstat {shlex.quote(str(external_id))}"
        try:
            result = self.run(command, timeout=timeout)
        except (OSError, subprocess.SubprocessError) as error:
            return {"state": "unknown", "detail": str(error)}
        text = ((result.stdout or "") + (result.stderr or "")).strip()
        return {"state": classify_state(text, scheduler),
                "detail": text[-400:]}

    def exists(self, path, timeout=30, dry_run=False):
        """Return True when ``path`` exists on this host."""
        if dry_run:
            return False
        if self.transport == "local":
            return os.path.exists(str(path))
        result = self.run(f"test -e {shlex.quote(str(path))}", timeout=timeout)
        return result.returncode == 0

    def start_background(self, command, cwd=None, log_path=None, dry_run=False):
        """Start ``command`` detached and return its pid (as a string).

        Used by the ``process`` execution mode for long runs that must survive
        the automag process (``--once``).
        """
        line = self.wrap(command, cwd=cwd)
        if log_path:
            line = f"{line} >> {shlex.quote(str(log_path))} 2>&1"
        if dry_run:
            print(f"    [dry-run @{self.name}] {line} &")
            return "dry-run"
        if self.transport == "local":
            # A new session keeps the child alive after automag exits.
            process = subprocess.Popen(
                ["/bin/bash", "-c", f"{line} & echo $!"],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                start_new_session=True)
            out, err = process.communicate(timeout=60)
            pid = (out or "").strip().splitlines()[-1:] or [""]
            return pid[0]
        result = self.run(f"{line} > /dev/null 2>&1 & echo $!", timeout=60)
        return ((result.stdout or "").strip().splitlines() or [""])[-1]

    def check(self, timeout=20, dry_run=False):
        """Verify that the host is reachable; return an error message or None."""
        if dry_run:
            return None
        try:
            result = self.run("echo automag-ok", timeout=timeout)
        except (OSError, subprocess.SubprocessError) as error:
            return str(error)
        text = (result.stdout or "").strip()
        if result.returncode != 0 or "automag-ok" not in text:
            detail = (result.stderr or text or "no output").strip()
            return detail.splitlines()[-1] if detail else "unreachable"
        return None


def classify_state(text, scheduler=None):
    """Map scheduler output to ``pending``/``running``/``finished``/``failed``.

    ``unknown`` is returned when the output carries no recognisable state (the
    job is no longer known to the scheduler).
    """
    if not text:
        return "unknown"
    lowered = text.lower()
    if scheduler == "slurm":
        for word, state in (
                ("completed", "finished"), ("running", "running"),
                ("completing", "running"), ("pending", "pending"),
                ("configuring", "pending"), ("requeued", "pending"),
                ("suspended", "pending"),
                ("failed", "failed"), ("cancelled", "failed"),
                ("timeout", "failed"), ("node_fail", "failed"),
                ("out_of_memory", "failed"), ("preempted", "failed")):
            if word in lowered:
                return state
        if "invalid job id" in lowered:
            return "unknown"
        return "unknown"
    if scheduler in ("pbs", "pbspro", "torque"):
        if "unknown job" in lowered or "unknown job id" in lowered:
            return "unknown"
        for token, state in ((" C ", "finished"), (" R ", "running"),
                             (" E ", "running"), (" Q ", "pending"),
                             (" H ", "pending"), (" W ", "pending"),
                             (" S ", "pending")):
            if token in f" {text} ":
                return state
        if "job_state = c" in lowered:
            return "finished"
        return "unknown"
    # process / pid: ``ps -p <pid> -o stat=`` prints nothing when it is gone
    return "running" if text.strip() else "finished"


def parse_job_id(text, scheduler=None):
    """Extract a scheduler job id from ``sbatch``/``qsub`` output or a pid."""
    import re

    if not text:
        return None
    for pattern in (r"Submitted batch job\s+(\d+)",   # slurm
                    r"^\s*(\d+)\[\S*\]?\s*$",        # pbs array
                    r"^\s*(\d+)\s*$"):               # pbs / pid
        match = re.search(pattern, text, re.MULTILINE)
        if match:
            return match.group(1)
    match = re.search(r"\b(\d{3,})\b", text)
    return match.group(1) if match else None


def _completed(returncode, stdout, stderr):
    return subprocess.CompletedProcess([], returncode, stdout, stderr)


def resolve_hosts(config, local_workdir=None):
    """Return ``{name: Host}`` merged from the config, always with ``local``."""
    from common.paths import code_root

    hosts = {"local": Host(name="local", transport="local")}
    for name, value in (config.hosts or {}).items():
        if name == "local":
            # Allow overriding the local host (paths, environment).
            hosts["local"] = Host.from_config("local", value)
            continue
        hosts[name] = Host.from_config(name, value)
    resolved = {name: host.resolved(local_workdir)
                for name, host in hosts.items()}
    for host in resolved.values():
        # A 'shell' host runs on this machine (inside a container, under
        # 'nice', ...), so the automag checkout is the local one unless the
        # user says otherwise.  An 'ssh' host needs an explicit path.
        if host.transport == "shell" and not host.automag_home and not host.automag_exe:
            host.automag_home = str(code_root())
    return resolved


def host_for(instance, hosts):
    """Return the host of one stage instance (named reference or inline)."""
    host = instance.config.get("host")
    if host is None:
        host_cfg = instance.host
    else:
        host_cfg = host
    if host_cfg is None:
        return hosts["local"]
    if isinstance(host_cfg, str):
        if host_cfg not in hosts:
            known = ", ".join(sorted(hosts))
            raise HostError(
                f"Instance '{instance.id}' references unknown host "
                f"'{host_cfg}'. Known hosts: {known}.")
        return hosts[host_cfg]
    inline = Host.from_config(f"{instance.id}:inline", host_cfg)
    if inline.workdir is None:
        inline.workdir = hosts["local"].workdir
    return inline


def hosts_summary(hosts):
    return "\n".join(host.describe() for host in hosts.values())
