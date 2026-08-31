"""Binding against activation for flagellin at TLR5."""

from __future__ import annotations

from typing import Any

from . import _core

__all__ = ["decompose"]


def decompose(observed_fold: float) -> dict[str, Any]:
    """Split an observed activation ratio into a 1:1 part and a 2:2 part, for
    the measured `StFliC`/`RhFlaB` construct pair.

    Returns a dict with `observed_fold`, `monomer_fold`, `dimer_fold` and
    `proofreading_steps`.
    """
    return _core.receptor_decompose(observed_fold)
