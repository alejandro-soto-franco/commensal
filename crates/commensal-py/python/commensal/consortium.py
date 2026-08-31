"""Colonisation and butyrate production for a defined gut consortium."""

from __future__ import annotations

from typing import Any

from numpy.typing import NDArray

from . import _core

__all__ = ["steady_state", "simulate"]


def steady_state(input: float) -> dict[str, Any]:
    """Where a constant substrate supply `input` leaves the lumen once
    transients have gone, for the reference *A. caccae* consortium.

    Returns a dict with `substrate`, `biomass`, `butyrate` and `persists`,
    the last saying whether the strain withstands washout.
    """
    return _core.consortium_steady_state(input)


def simulate(
    dose: float, hours: float, dt: float
) -> tuple[NDArray[Any], NDArray[Any]]:
    """Integrate the chemostat under a constant substrate supply `dose`, for
    `hours` at step `dt`, from a fixed inoculum.

    Returns `(times, states)`, where `states` has one row per time point and
    three columns: substrate, biomass, butyrate.

    Raises `ValueError` if `hours` or `dt` is not positive.
    """
    return _core.consortium_simulate(dose, hours, dt)
