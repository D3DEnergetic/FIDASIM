"""Helpers for generating regression-test reference distributions."""

from .config import read_config
from .workflow import generate_outputs

__all__ = ["generate_outputs", "read_config"]
