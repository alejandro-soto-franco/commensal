import pytest

import commensal
from commensal import amplicon, consortium, receptor


def test_dereplication_counts_and_orders():
    reads = [b"ACGT", b"ACGT", b"TTTT", b"ACGT", b"GGGG", b"TTTT"]
    got = amplicon.dereplicate(reads)
    assert got[0] == (b"ACGT", 3)
    assert got[1] == (b"TTTT", 2)
    assert got[2] == (b"GGGG", 1)


def test_bimera_flags_a_join_of_two_abundant_parents():
    uniques = [(b"AAAACCCC", 40), (b"GGGGTTTT", 40), (b"AAAATTTT", 2)]
    flags = amplicon.flag_bimeras(uniques)
    assert flags == [False, False, True]


def test_steady_state_persists_above_the_washout_threshold():
    below = consortium.steady_state(0.002)
    above = consortium.steady_state(0.5)
    assert below["persists"] is False
    assert above["persists"] is True
    assert above["butyrate"] > 0.0


def test_decomposition_splits_an_observed_fold():
    d = receptor.decompose(94.0)
    assert d["monomer_fold"] * d["dimer_fold"] == pytest.approx(94.0, rel=1e-9)


def test_version_is_exposed():
    assert commensal.__version__


def test_simulate_rejects_a_non_positive_step():
    with pytest.raises(ValueError):
        consortium.simulate(0.5, 0.0, 0.01)
    with pytest.raises(ValueError):
        consortium.simulate(0.5, 24.0, 0.0)
