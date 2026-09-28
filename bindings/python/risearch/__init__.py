"""Python interface to RIsearch."""

from ._api import index, search
from ._native import TargetRegistry

__all__ = [
    "TargetRegistry",
    "index",
    "search",
]
