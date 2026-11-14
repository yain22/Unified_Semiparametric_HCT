"""Shared output locations for the planning analyses."""

from pathlib import Path
import os


def output_root():
    default = Path(__file__).resolve().parents[2] / "results" / "planning"
    return Path(os.environ.get("HCT_PLANNING_OUTPUT", default)).resolve()


def output_dir(name):
    path = output_root() / name
    path.mkdir(parents=True, exist_ok=True)
    return path
