"""Type stubs for the compiled extension."""

from typing import Any

from numpy.typing import NDArray

__version__: str

def amplicon_dereplicate(reads: list[bytes]) -> list[tuple[bytes, int]]: ...
def amplicon_flag_bimeras(
    uniques: list[tuple[bytes, int]],
    min_fold: float = ...,
    min_parent: int = ...,
    max_shift: int = ...,
) -> list[bool]: ...
def consortium_steady_state(input: float) -> dict[str, Any]: ...
def consortium_simulate(
    dose: float, hours: float, dt: float
) -> tuple[NDArray[Any], NDArray[Any]]: ...
def receptor_decompose(observed_fold: float) -> dict[str, Any]: ...
