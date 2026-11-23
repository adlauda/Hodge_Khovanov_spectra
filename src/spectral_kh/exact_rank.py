"""Eigensolver-independent Betti numbers and rank-indexed spectral gaps.

Why this module exists
----------------------
An early database driver asked ARPACK for the 50 smallest-magnitude eigenvalues
(``eigsh(..., which="SM")``) and counted those below a threshold.  Lanczos does
not reliably return every copy of a repeated eigenvalue, and that solver
undercounted the kernel of the Hodge Laplacian in 1,774 stored bidegrees of
the knot database, while still reporting the correct smallest positive
eigenvalue.  See ``Paper/Current/reproducibility/NUMERICAL-REVIEW.md``.

The functions here never decide by a threshold whether an eigenvalue is zero.

* ``betti_mod_p`` returns ``dim - rank(d_out) - rank(d_in)`` with ranks taken
  modulo a prime, after checking ``d_out @ d_in == 0`` in integer arithmetic.
  A modular rank never exceeds the rational rank, so the result is a rigorous
  *upper* bound on the rational Betti number.  It equals the rational Betti
  number unless the prime divides the order of a torsion class of the
  integral homology in an adjacent degree.
* ``gram_floor`` reads the smallest positive eigenvalue of a differential's
  Gram matrix at the index fixed by the exact rank.  By the Hodge
  decomposition the positive spectrum of the Laplacian in degree h is the
  union, with multiplicity, of the positive spectra of ``d_out^T d_out`` and
  ``d_in d_in^T``.
* ``low_spectrum`` combines the two and is what
  ``KhovanovEngine.get_low_spectrum_adaptive`` calls.

Only NumPy and SciPy are required.  ``self_test`` additionally uses
python-flint to compare the modular rank with an independent implementation.
"""

from __future__ import annotations

import numpy as np
import scipy.linalg
from scipy import sparse
from scipy.sparse.linalg import eigsh

PRIME = 65521  # largest prime below 2**16; see rank_mod_p_dense

# Largest differential (rows or columns) handled with dense exact ranks.
# A dense float64 copy of a 25000 x 25000 matrix needs 5 GB.
RANK_MAX = 25000
# Largest Gram matrix diagonalised densely when reading the gap.
GRAM_MAX = 3000
# A value at an index the exact rank says is zero must be below this fraction
# of the gap for an iterative solve to be accepted.
LEAK_MAX = 1e-6


class BlockTooLarge(ValueError):
    """The bidegree's differentials exceed the dense exact-rank limit."""


def _reduce(x, prime):
    """Residues in [0, prime) of an array of exact integers held in float64.

    ``x / prime`` is correctly rounded, so for ``|x| < 2**52`` its floor is the
    true integer quotient: at a multiple of the prime the quotient is exact,
    and elsewhere it is at least ``1/prime`` from an integer, far more than
    the rounding error.
    """
    return x - prime * np.floor(x / prime)


def rank_mod_p_dense(a, prime=PRIME, block=64):
    """Exact rank over GF(prime) of an integer array, by blocked elimination.

    The array holds exact integers in float64.  Each panel of ``block`` columns
    is reduced modulo the prime and eliminated; the remaining columns are then
    updated with one matrix product.  Columns in which no remaining row is
    nonzero are skipped; they stay zero in the remaining rows, because later
    pivot rows of the same panel were themselves zero in that column.

    The trailing block is *not* reduced after each update.  Every update adds
    at most ``block * prime**2`` in absolute value, so after all panels the
    entries stay below ``(columns + block) * prime**2``, which the assertion
    keeps under 2**52.  Entries are reduced when they are read: as a panel, or
    as pivot rows before they multiply anything.
    """
    a = _reduce(np.asarray(a, dtype=np.float64), prime)
    if a.shape[0] < a.shape[1]:
        a = np.ascontiguousarray(a.T)
    m, n = a.shape
    assert (n + 2 * block) * prime * prime < 2 ** 52
    r = 0
    for j0 in range(0, n, block):
        if r == m:
            break
        j1 = min(j0 + block, n)
        width = j1 - j0
        panel = _reduce(a[r:, j0:j1], prime)
        rows = panel.shape[0]
        mult = np.zeros((rows, width))
        k = 0
        for c in range(width):
            if k == rows:
                break
            hit = np.flatnonzero(panel[k:, c])
            if hit.size == 0:
                continue
            piv = k + int(hit[0])
            if piv != k:
                panel[[k, piv]] = panel[[piv, k]]
                mult[[k, piv]] = mult[[piv, k]]
                a[[r + k, r + piv], j1:] = a[[r + piv, r + k], j1:]
            inv = pow(int(panel[k, c]), -1, prime)
            col = _reduce(panel[k + 1:, c] * inv, prime)
            mult[k + 1:, k] = col
            panel[k + 1:, c:] = _reduce(panel[k + 1:, c:] - np.outer(col, panel[k, c:]), prime)
            k += 1
        if k and j1 < n:
            low = mult[:, :k]
            top = _reduce(a[r:r + k, j1:], prime)
            for t in range(1, k):
                top[t] = _reduce(top[t] - low[t, :t] @ top[:t], prime)
            a[r:r + k, j1:] = top
            if rows > k:
                a[r + k:, j1:] -= low[k:] @ top
        r += k
    return r


def rank_mod_p(matrix, prime=PRIME):
    """Exact rank over GF(prime) of a scipy sparse integer matrix."""
    matrix = sparse.csr_matrix(matrix)
    if matrix.nnz == 0:
        return 0
    # Zero rows and columns do not affect the rank.
    matrix = matrix[np.diff(matrix.indptr) > 0]
    csc = matrix.tocsc()
    csc = csc[:, np.diff(csc.indptr) > 0]
    return rank_mod_p_dense(csc.toarray(), prime)


def rank_flint(a, prime=PRIME):
    """Reference implementation used only by the self-test."""
    import flint
    a = np.mod(np.asarray(a, dtype=np.int64), prime)
    return flint.nmod_mat(a.shape[0], a.shape[1], a.ravel().tolist(), prime).rank()


def self_test(trials=60, seed=0):
    """Compare the blocked elimination with FLINT on rank-deficient matrices."""
    rng = np.random.default_rng(seed)
    for trial in range(trials):
        m, n = rng.integers(1, 90, size=2)
        true_rank = int(rng.integers(0, min(m, n) + 1))
        if trial % 3 == 0:  # sparse +-1 entries, like a Khovanov differential
            a = rng.choice([-1, 0, 0, 0, 0, 1], size=(m, n))
        else:  # product of two factors, so the rank is at most true_rank
            a = rng.integers(-3, 4, size=(m, true_rank)) @ rng.integers(-3, 4, size=(true_rank, n))
        for block in (1, 7, 64):
            got = rank_mod_p_dense(a, PRIME, block=block)
            want = rank_flint(a)
            assert got == want, (trial, m, n, block, got, want)
    print(f"self-test passed: {trials} matrices, three block sizes, agree with FLINT")


def gram_floor(matrix, rank):
    """Smallest positive eigenvalue of matrix^T matrix, indexed by exact rank.

    Returns ``None`` for rank 0.  ``floor`` is the eigenvalue of the smaller
    Gram matrix at ascending index ``k - rank``, ``next`` the one above it,
    and ``zero_leak`` the magnitude of the one below it, which is zero in
    exact arithmetic and so measures the floating-point error.
    """
    if rank == 0:
        return None
    m, n = matrix.shape
    gram = (matrix @ matrix.T) if m <= n else (matrix.T @ matrix)
    k = gram.shape[0]
    nullity = k - rank
    lo, hi = max(nullity - 1, 0), min(nullity + 1, k - 1)
    w = scipy.linalg.eigvalsh(gram.toarray().astype(float), subset_by_index=[lo, hi])
    values = dict(zip(range(lo, hi + 1), w))
    return {
        "gram_dim": k,
        "floor": float(values[nullity]),
        "next": float(values[nullity + 1]) if nullity + 1 in values else None,
        "zero_leak": float(abs(values[nullity - 1])) if nullity >= 1 else 0.0,
    }


def differentials(engine, h, q):
    """Outgoing and incoming integer differentials of bidegree (h, q)."""
    d_out = engine.get_differential_matrix(h, q).astype(np.int64)
    d_in = engine.get_differential_matrix(h - 1, q).astype(np.int64)
    return d_out, d_in


def is_exactly_zero(product):
    return bool(product.nnz == 0 or not product.data.any())


def betti_mod_p(engine, h, q, prime=PRIME, rank_max=RANK_MAX):
    """Exact-rank Betti number of bidegree (h, q); see the module docstring.

    Returns ``(betti, rank_out, rank_in)``.  Raises ``ArithmeticError`` if the
    differentials do not compose to zero, and ``BlockTooLarge`` if either one
    has more than ``rank_max`` rows or columns.
    """
    dim = len(engine.chain_groups.get((h, q), []))
    d_out, d_in = differentials(engine, h, q)
    if max(d_out.shape + d_in.shape) > rank_max:
        raise BlockTooLarge(f"bidegree ({h}, {q}): differential shapes "
                            f"{d_out.shape} and {d_in.shape} exceed {rank_max}")
    if not is_exactly_zero(d_out @ d_in):
        raise ArithmeticError(f"d^2 != 0 in bidegree ({h}, {q})")
    r_out, r_in = rank_mod_p(d_out, prime), rank_mod_p(d_in, prime)
    return dim - r_out - r_in, r_out, r_in


def _laplacian_values(engine, h, q, betti, count):
    """Eigenvalues ``betti .. betti+count-1`` of the Laplacian by shift-invert.

    The solver is asked for ``betti + count + 1`` eigenvalues nearest a small
    negative shift.  Had it missed a copy of the zero eigenvalue, the value at
    index ``betti - 1`` would be a positive eigenvalue instead of rounding
    noise, so the result is accepted only if that value is below ``LEAK_MAX``
    times the gap.  Returns ``None`` when the check fails.
    """
    lap = engine.compute_laplacian(h, q).astype(float)
    dim = lap.shape[0]
    k = betti + count + 1
    if k >= dim - 1:
        w = scipy.linalg.eigvalsh(lap.toarray())
    else:
        w = np.sort(eigsh(lap, k=k, sigma=-1e-5, return_eigenvectors=False))
        if betti and abs(w[betti - 1]) > LEAK_MAX * w[betti]:
            return None
    return [float(x) for x in w[betti:betti + count]]


def low_spectrum(engine, h, q, count=2, prime=PRIME, rank_max=RANK_MAX, gram_max=GRAM_MAX):
    """Betti number and the ``count`` smallest positive Laplacian eigenvalues.

    The Betti number comes from exact modular ranks.  The eigenvalues are read
    at the index the ranks dictate: from dense Gram matrices of the two
    differentials when those are at most ``gram_max`` on a side, otherwise
    from a shift-invert solve of the Laplacian that is checked against the
    exact Betti number and replaced by the dense computation if the check
    fails.  ``count`` may be 1 or 2.
    """
    if count not in (1, 2):
        raise ValueError("count must be 1 or 2")
    dim = len(engine.chain_groups.get((h, q), []))
    if dim == 0:
        return 0, []
    betti, r_out, r_in = betti_mod_p(engine, h, q, prime, rank_max)
    if betti == dim:
        return betti, []
    d_out, d_in = differentials(engine, h, q)
    gram_dims = [min(d.shape) for d, r in ((d_out, r_out), (d_in, r_in)) if r]
    if max(gram_dims) > gram_max:
        values = _laplacian_values(engine, h, q, betti, count)
        if values is not None:
            return betti, values
    values = []
    for d, r in ((d_out, r_out), (d_in, r_in)):
        f = gram_floor(d, r)
        if f is not None:
            values += [f["floor"]] + ([f["next"]] if f["next"] is not None else [])
    return betti, sorted(values)[:count]


if __name__ == "__main__":
    self_test()
