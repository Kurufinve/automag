"""
automag.pipeline.cli
====================

Command line entry point of the unified pipeline.

Run modes:

``(default)``
    execute the requested phases of every instance once, in order.
``--wait``
    run the whole pipeline continuously, waiting for the calculations between
    the steps (state kept in SQLite so an interrupted run can be resumed).
``--once``
    advance the pipeline until a step has to wait for external programs,
    record everything and exit.  Run the same command again to continue.
``--daemon``
    like ``--wait`` but never exits while there is pending work.
``--status`` / ``--report``
    print the stored state and exit.

.. codeauthor:: Automag contributors
"""
from __future__ import annotations


import argparse
import subprocess
import sys

from pipeline.config import ConfigError, PipelineConfig
from pipeline.runner import PipelineRunner, stage_summary


def build_parser():
    parser = argparse.ArgumentParser(
        prog="run_automag",
        description="Run the automag pipeline from a single input_automag.py file.",
    )
    parser.add_argument(
        "config",
        nargs="?",
        default=None,
        help=("controlling file or directory (default: input_automag.py, "
              "then input_automag.yaml, then input_automag.json in the "
              "current directory)"),
    )
    parser.add_argument(
        "--steps",
        default=None,
        help="override the phases to run, e.g. 'submit' or 'submit,postprocess'",
    )
    parser.add_argument(
        "--stages",
        default=None,
        help="run only the given stages, e.g. 'coll,relax' (comma separated)",
    )
    parser.add_argument(
        "--instances",
        default=None,
        help=("run only the given stage instance ids, e.g. "
              "'coll_prophet,mae_vasp' (comma separated)"),
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--wait",
        action="store_true",
        help=("run the complete pipeline and wait for submitted calculations "
              "between phases, so no manual step is needed (state is stored "
              "in SQLite and the run can be resumed)"),
    )
    mode.add_argument(
        "--once",
        action="store_true",
        help=("advance the pipeline until a step has to wait for external "
              "programs, record the started jobs and exit; rerun to continue"),
    )
    mode.add_argument(
        "--daemon",
        action="store_true",
        help="like --wait but keep running to supervise the started jobs",
    )
    mode.add_argument(
        "--status",
        action="store_true",
        help="print the stored pipeline state and exit",
    )
    parser.add_argument(
        "--report",
        action="store_true",
        help="print a detailed state report (steps and running jobs) and exit",
    )
    parser.add_argument(
        "--execution",
        default=None,
        choices=("auto", "inprocess", "process", "scheduler"),
        help=("force the execution mode of every step: auto (default), "
              "inprocess, process or scheduler"),
    )
    parser.add_argument(
        "--fresh",
        action="store_true",
        help="ignore / delete the stored state and start the pipeline over",
    )
    parser.add_argument(
        "--state-file",
        default=None,
        help="path of the SQLite state file (default: <workdir>/automag_state.sqlite)",
    )
    parser.add_argument(
        "--poll-interval",
        type=float,
        default=None,
        help="polling interval in seconds for --wait/--once (default from input file)",
    )
    parser.add_argument(
        "--wait-timeout",
        type=float,
        default=None,
        help="maximum time in seconds to wait for one phase (default from input file)",
    )
    parser.add_argument(
        "--wait-interval",
        type=float,
        default=None,
        help="polling interval in seconds for --wait (default from input file)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="prepare directories and input.py files but do not execute anything",
    )
    parser.add_argument(
        "--list-stages",
        action="store_true",
        help="print the available stages and phases and exit",
    )
    parser.add_argument(
        "--list-hosts",
        action="store_true",
        help="print the resolved execution hosts and exit",
    )
    return parser


def _warn_unavailable_codes(config):
    """Point out external codes the build step could not provide."""
    from common import external

    overrides = config.common.get("external_codes") or None
    broken = [name for name in external.CODES
              if not external.available(name, overrides)]
    if not broken:
        return
    print("WARNING: these external codes are not available: "
          + ", ".join(broken))
    print("         The pipeline steps that need them will stop with an error.")
    print("         Fix with './build.sh --external' or by passing the paths, "
          "e.g. './build.sh --external --espins-path /path/to/mc.x'.")


def _print_state(config, state_file=None, report=False):
    from pipeline.state import StateStore
    from pipeline.runner import PipelineRunner

    runner = PipelineRunner(config, dry_run=True, mode="steps",
                            state_file=state_file)
    path = runner.state_path()
    if not path.is_file():
        print(f"No stored state at {path}.")
        return 1
    with StateStore(path, create=False) as state:
        print(state.report() if report else state.summary())
    return 0


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.list_stages:
        print(stage_summary())
        return 0

    try:
        config = PipelineConfig.load(args.config)
        if args.list_hosts:
            from pipeline import hosts as hosts_module

            resolved = hosts_module.resolve_hosts(
                config, local_workdir=config.workdir)
            print(hosts_module.hosts_summary(resolved))
            return 0
        if args.status or args.report:
            return _print_state(config, args.state_file, report=args.report)
        _warn_unavailable_codes(config)
        stages = None
        if args.stages:
            stages = [item.strip() for item in args.stages.split(",") if item.strip()]
        instances = None
        if args.instances:
            instances = [item.strip() for item in args.instances.split(",")
                         if item.strip()]
        mode = "steps"
        if args.wait:
            mode = "wait"
        elif args.once:
            mode = "once"
        elif args.daemon:
            mode = "daemon"
        runner = PipelineRunner(
            config, steps=args.steps, stages=stages, instances=instances,
            dry_run=args.dry_run, mode=mode, execution=args.execution,
            fresh=args.fresh, state_file=args.state_file,
            poll_interval=args.poll_interval,
            wait_timeout=args.wait_timeout,
            wait_interval=args.wait_interval,
        )
        runner.run()
    except ConfigError as error:
        print(f"Configuration error: {error}", file=sys.stderr)
        return 2
    except subprocess.CalledProcessError as error:
        print(f"Stage failed: {error}", file=sys.stderr)
        return 1

    print("\nPipeline finished.")
    return 0
