"""Reads to unique sequences."""

from __future__ import annotations

from . import _core

__all__ = ["dereplicate", "flag_bimeras"]


def dereplicate(reads: list[bytes]) -> list[tuple[bytes, int]]:
    """Collapse reads to unique sequences with counts, most abundant first.

    Ties break on the sequence itself, so a run repeats.
    """
    return _core.amplicon_dereplicate(reads)


def flag_bimeras(
    uniques: list[tuple[bytes, int]],
    min_fold: float = 2.0,
    min_parent: int = 8,
    max_shift: int = 16,
) -> list[bool]:
    """Flag every sequence that reads as the left of one parent and the right
    of another. The defaults are the reference implementation's."""
    return _core.amplicon_flag_bimeras(uniques, min_fold, min_parent, max_shift)
