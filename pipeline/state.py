"""
automag.pipeline.state
======================

Persistent pipeline state (SQLite) so that automag can be interrupted and
resumed, and so that individual calculations can be retried selectively.

Three levels are recorded:

``steps``
    one row per ``(instance id, phase)`` - the unit the runner schedules.
    Status: ``pending`` / ``running`` / ``waiting`` / ``done`` / ``failed`` /
    ``skipped``.  ``config_hash`` stores the hash of the settings the step ran
    with, which is how a changed input file is detected.

``tasks``
    one row per *submitted external job*: a scheduler job id, a PID or a remote
    pid, plus the list of calculation items it was supposed to produce.

``items``
    one row per individual calculation (one magnetic configuration).  This is
    what makes a partial restart possible: after a batch job dies, only the
    items without a result are resubmitted, and every item counts its own
    attempts.

The database lives next to the results (``<workdir>/automag_state.sqlite`` by
default) and is written with ``sqlite3`` from the standard library, so it
survives a killed automag process and needs no extra service.

.. codeauthor:: Automag contributors
"""
from __future__ import annotations


import json
import sqlite3
import time
import uuid
from contextlib import contextmanager
from pathlib import Path

SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (
    key   TEXT PRIMARY KEY,
    value TEXT
);

CREATE TABLE IF NOT EXISTS steps (
    step_key    TEXT PRIMARY KEY,
    pipeline    TEXT,
    instance    TEXT,
    stage       TEXT,
    phase       TEXT,
    status      TEXT NOT NULL DEFAULT 'pending',
    attempts    INTEGER NOT NULL DEFAULT 0,
    host        TEXT,
    config_hash TEXT,
    started_at  REAL,
    updated_at  REAL,
    finished_at REAL,
    reset_at    REAL,
    detail      TEXT
);

CREATE TABLE IF NOT EXISTS tasks (
    task_key    TEXT PRIMARY KEY,
    step_key    TEXT NOT NULL,
    host        TEXT,
    kind        TEXT,
    external_id TEXT,
    workdir     TEXT,
    command     TEXT,
    status      TEXT NOT NULL DEFAULT 'pending',
    attempts    INTEGER NOT NULL DEFAULT 0,
    items       TEXT,
    launcher    TEXT,
    updated_at  REAL,
    detail      TEXT
);

CREATE TABLE IF NOT EXISTS items (
    item_key TEXT PRIMARY KEY,
    step_key TEXT NOT NULL,
    path     TEXT,
    status   TEXT NOT NULL DEFAULT 'pending',
    attempts INTEGER NOT NULL DEFAULT 0,
    task_key TEXT,
    updated_at REAL,
    detail   TEXT
);

CREATE INDEX IF NOT EXISTS idx_items_step ON items(step_key, status);
CREATE INDEX IF NOT EXISTS idx_tasks_step ON tasks(step_key, status);
"""

ACTIVE_STATUSES = ("pending", "running", "waiting", "submitted")

#: Columns added after the first release; ``CREATE TABLE IF NOT EXISTS`` does
#: not touch an existing database, so they are migrated explicitly.
MIGRATIONS = (
    ("steps", "reset_at", "REAL"),
    ("tasks", "launcher", "TEXT"),
)


def now():
    return time.time()


class StateStore:
    """Thin, explicit wrapper around the sqlite database."""

    def __init__(self, path, create=True):
        self.path = Path(path).expanduser()
        self.create = create
        self._connection = None

    # ------------------------------------------------------------------
    # Connection
    # ------------------------------------------------------------------
    @property
    def connection(self):
        if self._connection is None:
            if self.create:
                self.path.parent.mkdir(parents=True, exist_ok=True)
            if not self.path.is_file() and not self.create:
                raise FileNotFoundError(
                    f"No automag state file at {self.path}. Run the pipeline "
                    "once, or pass '--fresh' to ignore the stored state.")
            self._connection = sqlite3.connect(str(self.path))
            self._connection.row_factory = sqlite3.Row
            self._connection.executescript(SCHEMA)
            self._migrate()
            self._connection.commit()
        return self._connection

    def _migrate(self):
        """Add columns that older databases do not have yet."""
        for table, column, kind in MIGRATIONS:
            existing = {row["name"] for row in self._connection.execute(
                f"PRAGMA table_info({table})")}
            if column not in existing:
                self._connection.execute(
                    f"ALTER TABLE {table} ADD COLUMN {column} {kind}")

    def close(self):
        if self._connection is not None:
            self._connection.commit()
            self._connection.close()
            self._connection = None

    def __enter__(self):
        return self

    def __exit__(self, *exception):
        self.close()

    @contextmanager
    def transaction(self):
        connection = self.connection
        try:
            yield connection
            connection.commit()
        except Exception:
            connection.rollback()
            raise

    # ------------------------------------------------------------------
    # meta
    # ------------------------------------------------------------------
    def get_meta(self, key, default=None):
        row = self.connection.execute(
            "SELECT value FROM meta WHERE key = ?", (key,)).fetchone()
        return row["value"] if row else default

    def set_meta(self, key, value):
        with self.transaction() as connection:
            connection.execute(
                "INSERT INTO meta(key, value) VALUES(?, ?) "
                "ON CONFLICT(key) DO UPDATE SET value = excluded.value",
                (key, None if value is None else str(value)))

    # ------------------------------------------------------------------
    # steps
    # ------------------------------------------------------------------
    def get_step(self, step_key):
        row = self.connection.execute(
            "SELECT * FROM steps WHERE step_key = ?", (step_key,)).fetchone()
        return dict(row) if row else None

    def step_status(self, step_key):
        step = self.get_step(step_key)
        return step["status"] if step else None

    def steps(self):
        rows = self.connection.execute("SELECT * FROM steps").fetchall()
        return {row["step_key"]: dict(row) for row in rows}

    def upsert_step(self, step_key, **fields):
        fields.setdefault("updated_at", now())
        detail = fields.get("detail")
        if isinstance(detail, (dict, list)):
            fields["detail"] = json.dumps(detail, default=str)
        existing = self.get_step(step_key)
        with self.transaction() as connection:
            if existing is None:
                connection.execute(
                    "INSERT INTO steps(step_key, status, updated_at) "
                    "VALUES(?, ?, ?)",
                    (step_key, fields.pop("status", "pending"), now()))
                existing = {"step_key": step_key}
            assignments = ", ".join(f"{key} = ?" for key in fields)
            connection.execute(
                f"UPDATE steps SET {assignments} WHERE step_key = ?",
                (*fields.values(), step_key))
        return self.get_step(step_key)

    def mark_step(self, step_key, status, **fields):
        fields["status"] = status
        if status in ("done", "failed", "skipped"):
            fields.setdefault("finished_at", now())
        if status == "running":
            fields.setdefault("started_at", now())
        return self.upsert_step(step_key, **fields)

    def reset_steps(self, predicate=None):
        """Forget the stored state of the matching steps (and their children)."""
        for step_key, step in self.steps().items():
            if predicate is not None and not predicate(step_key, step):
                continue
            self.reset_step(step_key)

    def reset_step(self, step_key, reset_at=None, keep_record=False):
        """Forget one step: its tasks and items are deleted.

        With ``keep_record`` the step row itself survives with status
        ``pending`` and ``reset_at`` set, so that results produced before the
        reset are not mistaken for results of the new run.
        """
        with self.transaction() as connection:
            connection.execute("DELETE FROM items WHERE step_key = ?",
                               (step_key,))
            connection.execute("DELETE FROM tasks WHERE step_key = ?",
                               (step_key,))
            if keep_record:
                connection.execute(
                    "UPDATE steps SET status = 'pending', attempts = 0, "
                    "started_at = NULL, finished_at = NULL, reset_at = ?, "
                    "updated_at = ?, detail = NULL WHERE step_key = ?",
                    (now() if reset_at is None else reset_at, now(), step_key))
            else:
                connection.execute("DELETE FROM steps WHERE step_key = ?",
                                   (step_key,))

    # ------------------------------------------------------------------
    # items (individual calculations)
    # ------------------------------------------------------------------
    def add_items(self, step_key, paths, prefix=""):
        """Register calculation directories; existing rows are left alone."""
        added = []
        with self.transaction() as connection:
            for path in paths:
                item_key = f"{prefix}{path}"
                cursor = connection.execute(
                    "INSERT INTO items(item_key, step_key, path, status, "
                    "updated_at) VALUES(?, ?, ?, 'pending', ?) "
                    "ON CONFLICT(item_key) DO NOTHING",
                    (item_key, step_key, str(path), now()))
                if cursor.rowcount:
                    added.append(item_key)
        return added

    def add_item(self, step_key, item_key, path=None, status="pending"):
        """Insert (or keep) one item explicitly."""
        with self.transaction() as connection:
            cursor = connection.execute(
                "INSERT INTO items(item_key, step_key, path, status, "
                "updated_at) VALUES(?, ?, ?, ?, ?) "
                "ON CONFLICT(item_key) DO NOTHING",
                (item_key, step_key, str(path or item_key), status, now()))
            return bool(cursor.rowcount)

    def items(self, step_key, status=None):
        query = "SELECT * FROM items WHERE step_key = ?"
        parameters = [step_key]
        if status is not None:
            statuses = (status,) if isinstance(status, str) else tuple(status)
            query += " AND status IN (%s)" % ", ".join("?" * len(statuses))
            parameters.extend(statuses)
        rows = self.connection.execute(query, parameters).fetchall()
        return [dict(row) for row in rows]

    def pending_items(self, step_key, max_attempts=None):
        """Items that are not done yet and still have retries left."""
        result = []
        for item in self.items(step_key):
            if item["status"] == "done":
                continue
            if max_attempts is not None and item["attempts"] >= max_attempts:
                continue
            result.append(item)
        return result

    def set_item(self, item_key, **fields):
        fields.setdefault("updated_at", now())
        assignments = ", ".join(f"{key} = ?" for key in fields)
        with self.transaction() as connection:
            connection.execute(
                f"UPDATE items SET {assignments} WHERE item_key = ?",
                (*fields.values(), item_key))
        return self.get_item(item_key)

    def get_item(self, item_key):
        row = self.connection.execute(
            "SELECT * FROM items WHERE item_key = ?", (item_key,)).fetchone()
        return dict(row) if row else None

    def bump_item_attempts(self, item_keys, task_key=None):
        with self.transaction() as connection:
            for item_key in item_keys:
                connection.execute(
                    "UPDATE items SET attempts = attempts + 1, status = 'running',"
                    " task_key = ?, updated_at = ? WHERE item_key = ?",
                    (task_key, now(), item_key))

    def items_progress(self, step_key):
        items = self.items(step_key)
        done = sum(1 for item in items if item["status"] == "done")
        return done, len(items)

    # ------------------------------------------------------------------
    # tasks (submitted external jobs)
    # ------------------------------------------------------------------
    def add_task(self, step_key, host=None, kind=None, workdir=None,
                 command=None, items=None, status="pending", detail=None,
                 launcher=None, external_id=None):
        task_key = uuid.uuid4().hex[:12]
        with self.transaction() as connection:
            connection.execute(
                "INSERT INTO tasks(task_key, step_key, host, kind, workdir, "
                "command, status, attempts, items, launcher, external_id, "
                "updated_at, detail) "
                "VALUES(?, ?, ?, ?, ?, ?, ?, 1, ?, ?, ?, ?, ?)",
                (task_key, step_key, host, kind, workdir, command, status,
                 json.dumps(list(items or [])), launcher, external_id, now(),
                 None if detail is None else json.dumps(detail, default=str)))
        return task_key

    def update_task(self, task_key, **fields):
        fields.setdefault("updated_at", now())
        if isinstance(fields.get("detail"), (dict, list)):
            fields["detail"] = json.dumps(fields["detail"], default=str)
        if isinstance(fields.get("items"), (list, tuple)):
            fields["items"] = json.dumps(list(fields["items"]))
        assignments = ", ".join(f"{key} = ?" for key in fields)
        with self.transaction() as connection:
            connection.execute(
                f"UPDATE tasks SET {assignments} WHERE task_key = ?",
                (*fields.values(), task_key))
        return self.get_task(task_key)

    def get_task(self, task_key):
        row = self.connection.execute(
            "SELECT * FROM tasks WHERE task_key = ?", (task_key,)).fetchone()
        return self._task_row(row)

    def tasks(self, step_key=None, status=None):
        query = "SELECT * FROM tasks"
        clauses, parameters = [], []
        if step_key is not None:
            clauses.append("step_key = ?")
            parameters.append(step_key)
        if status is not None:
            statuses = (status,) if isinstance(status, str) else tuple(status)
            clauses.append("status IN (%s)" % ", ".join("?" * len(statuses)))
            parameters.extend(statuses)
        if clauses:
            query += " WHERE " + " AND ".join(clauses)
        rows = self.connection.execute(query, parameters).fetchall()
        return [self._task_row(row) for row in rows]

    @staticmethod
    def _task_row(row):
        if row is None:
            return None
        task = dict(row)
        try:
            task["items"] = json.loads(task["items"] or "[]")
        except ValueError:
            task["items"] = []
        return task

    # ------------------------------------------------------------------
    # Reporting
    # ------------------------------------------------------------------
    def summary(self):
        """Return a compact multi-line report of the stored state."""
        lines = []
        for step_key, step in sorted(self.steps().items()):
            done, total = self.items_progress(step_key)
            extra = f" items {done}/{total}" if total else ""
            host = f" @{step['host']}" if step.get("host") else ""
            lines.append(f"{step_key:32s} {step['status']:9s}{host}{extra}")
            for task in self.tasks(step_key,
                                   status=("pending", "submitted", "running")):
                lines.append(
                    f"    task {task['task_key']} {task['status']:9s} "
                    f"{task['kind'] or '-':8s} id={task['external_id'] or '-'}"
                    f" @{task['host'] or '-'}")
        return "\n".join(lines)

    def report(self):
        """Detailed report: steps, running tasks and per-step item counters."""
        lines = [self.summary() or "(no state recorded)"]
        if self.get_meta("config_path"):
            lines.append(f"config: {self.get_meta('config_path')}")
        if self.get_meta("updated_at"):
            lines.append(
                "updated: "
                + time.strftime("%Y-%m-%d %H:%M:%S",
                                time.localtime(float(self.get_meta("updated_at")))))
        return "\n".join(lines)

    def is_finished(self, step_key):
        return self.step_status(step_key) in ("done", "skipped")
