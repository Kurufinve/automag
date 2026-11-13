"""
automag.pipeline
================

Unified, single-command pipeline.

.. codeauthor:: Automag contributors
"""

from pipeline.config import ConfigError, PipelineConfig
from pipeline.runner import PipelineRunner, stage_summary
from pipeline.stages import STAGE_ORDER, STAGES, Stage

__all__ = [
    "ConfigError",
    "PipelineConfig",
    "PipelineRunner",
    "STAGES",
    "STAGE_ORDER",
    "Stage",
    "stage_summary",
]
