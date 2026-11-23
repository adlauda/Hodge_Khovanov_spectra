"""Regression tests for exact-rank Betti numbers and rank-indexed gaps.

Run with ``pytest`` or directly: ``python tests/test_exact_rank.py``.

The expected Betti tables are the rational Khovanov homology in the grading
used by ``databases/knot_research.db``; they agree with KnotInfo.  PD codes
are those of ``data/knotinfo.csv``.

Set ``SPECTRAL_KH_SLOW=1`` to include a 4000-dimensional bidegree of
``12n_803`` whose Betti number the original truncated eigensolver recorded as
5 instead of 7 (about ten seconds on eight cores).
"""

import ast
import csv
import os
import sys
import warnings
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from spectral_kh import KhovanovEngine, exact_rank  # noqa: E402

KNOTS = {
    "3_1": (
        [[1, 5, 2, 4], [3, 1, 4, 6], [5, 3, 6, 2]],
        {(0, 1): 1, (0, 3): 1, (2, 5): 1, (3, 9): 1},
        1.0,
    ),
    "4_1": (
        [[4, 2, 5, 1], [8, 6, 1, 5], [6, 3, 7, 4], [2, 7, 3, 8]],
        {(-2, -5): 1, (-1, -1): 1, (0, -1): 1, (0, 1): 1, (1, 1): 1, (2, 5): 1},
        1.0,
    ),
    "8_19": (  # non-alternating, the torus knot T(3,4)
        [[2, 14, 3, 13], [5, 11, 6, 10], [7, 15, 8, 14], [9, 5, 10, 4],
         [11, 7, 12, 6], [12, 2, 13, 1], [15, 9, 16, 8], [16, 4, 1, 3]],
        {(0, 5): 1, (0, 7): 1, (2, 9): 1, (3, 13): 1, (4, 11): 1, (4, 13): 1,
         (5, 15): 1, (5, 17): 1},
        0.200110824773,
    ),
    "9_42": (  # non-alternating, homologically thick
        [[1, 5, 2, 4], [5, 11, 6, 10], [3, 8, 4, 9], [9, 2, 10, 3], [16, 11, 17, 12],
         [14, 8, 15, 7], [6, 16, 7, 15], [18, 13, 1, 14], [12, 17, 13, 18]],
        {(-2, -7): 1, (-1, -3): 1, (0, -3): 1, (0, -1): 1, (0, 1): 1, (1, -1): 1,
         (1, 1): 1, (2, 3): 1, (3, 3): 1, (4, 7): 1},
        0.043242811485,
    ),
}


def _check_knot(name):
    pd_code, expected, expected_min_gap = KNOTS[name]
    engine = KhovanovEngine(pd_code)
    betti_table, min_gap = {}, np.inf
    for (h, q), basis in sorted(engine.chain_groups.items()):
        if not basis:
            continue
        betti, gaps = engine.get_low_spectrum_adaptive(h, q)
        # Same block through the iterative path, and through a dense
        # diagonalisation of the Laplacian read at the exact-rank index.
        betti_iter, gaps_iter = exact_rank.low_spectrum(engine, h, q, gram_max=0)
        dense = np.linalg.eigvalsh(engine.compute_laplacian(h, q).toarray())
        assert betti == betti_iter
        assert np.allclose(gaps, gaps_iter, rtol=1e-9, atol=0)
        assert np.allclose(gaps, dense[betti:betti + 2], rtol=1e-9, atol=0)
        assert np.all(np.abs(dense[:betti]) < 1e-9)
        if betti:
            betti_table[(h, q)] = int(betti)
        if len(gaps):
            min_gap = min(min_gap, gaps[0])
    assert betti_table == expected, name
    assert abs(min_gap - expected_min_gap) < 1e-10, (name, min_gap)


def test_trefoil():
    _check_knot("3_1")


def test_figure_eight():
    _check_knot("4_1")


def test_8_19():
    _check_knot("8_19")


def test_9_42():
    _check_knot("9_42")


def test_rank_mod_p_matches_rational_rank():
    rng = np.random.default_rng(0)
    for _ in range(40):
        m, n = rng.integers(1, 60, size=2)
        r = int(rng.integers(0, min(m, n) + 1))
        a = rng.integers(-3, 4, size=(m, r)) @ rng.integers(-3, 4, size=(r, n))
        assert exact_rank.rank_mod_p_dense(a) == np.linalg.matrix_rank(a)
        for block in (1, 7):
            assert exact_rank.rank_mod_p_dense(a, block=block) == np.linalg.matrix_rank(a)


def test_large_block_falls_back_with_warning():
    engine = KhovanovEngine(KNOTS["4_1"][0])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        betti, gaps = engine.get_low_spectrum_adaptive(0, -1, rank_max=1)
    assert betti == 1 and len(gaps) == 2
    assert any(issubclass(w.category, RuntimeWarning) for w in caught)


def test_khovanov_polynomial_string():
    engine = KhovanovEngine(KNOTS["3_1"][0])
    assert engine.compute_khovanov_polynomial() == "q^{1} + q^{3} + t^{2}q^{5} + t^{3}q^{9}"


def test_block_the_truncated_solver_undercounted():
    if not os.environ.get("SPECTRAL_KH_SLOW"):
        return
    with open(REPO / "data" / "knotinfo.csv", newline="") as f:
        pd_code = next(ast.literal_eval(r["PD Notation"].replace(";", ","))
                       for r in csv.DictReader(f) if r["Name"] == "12n_803")
    engine = KhovanovEngine(pd_code)
    assert len(engine.chain_groups[(4, 11)]) == 4000
    betti, gaps = engine.get_low_spectrum_adaptive(4, 11)
    assert betti == 7
    assert abs(gaps[0] - EXPECTED_12N_803_GAP) < 1e-9 * EXPECTED_12N_803_GAP


EXPECTED_12N_803_GAP = 0.3090431086896931  # stored value, unchanged by the repair


if __name__ == "__main__":
    tests = [(k, v) for k, v in sorted(globals().items()) if k.startswith("test_")]
    for test_name, test in tests:
        test()
        print("ok", test_name)
