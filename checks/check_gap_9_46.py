"""Certify the common gap in two bidegrees of the KnotInfo diagram of 9_46.

Run alongside independent_kh.py with Python 3.10+ and assertions enabled.
Sylvester's law of inertia and fraction-free elimination give exact counts
of eigenvalues below the two rational endpoints. Only integers are used
in the inertia calculation. The PD code below is the 9_46 row from the
knot_research.db input archived by adlauda/Hodge_Khovanov_spectra at commit
ccfce9a3eaa33f6663ff1e21b158d73aa34b3bea.
"""
from collections import defaultdict
from fractions import Fraction
from independent_kh import Diagram, Khovanov

PD = (
    (2, 10, 3, 9), (3, 14, 4, 15), (6, 17, 7, 18),
    (8, 11, 9, 12), (10, 2, 11, 1), (13, 4, 14, 5),
    (15, 13, 16, 12), (16, 7, 17, 8), (18, 5, 1, 6),
)
LOW = Fraction("0.0289984")
HIGH = Fraction("0.0289985")


def knot_complex():
    crossings = []
    ends = defaultdict(list)
    for i, cross in enumerate(PD):
        ports = [(i, j) for j in range(4)]
        # The database's 0-smoothing pairs the first two and last two PD ports.
        crossings.append(tuple(ports[1:] + ports[:1]))
        for arc, port in zip(cross, ports):
            ends[arc].append(port)
    return Khovanov(Diagram(crossings, list(ends.values()), npos=3, nneg=6))


def laplacian(incoming, outgoing):
    dim = incoming.nrows
    matrix = [[0] * dim for _ in range(dim)]
    rows = defaultdict(dict)
    for i, column in enumerate(outgoing.columns):
        for j, coefficient in column.items():
            rows[j][i] = coefficient
    # d_in d_in^T sums outer products of columns, d_out^T d_out of rows.
    for v in incoming.columns + list(rows.values()):
        for i, x in v.items():
            for j, y in v.items():
                assert x.denominator == y.denominator == 1
                matrix[i][j] += int(x * y)
    return matrix


def count_below(matrix, threshold):
    """Count negative eigenvalues of b*matrix-a*I at threshold a/b."""
    a, b = threshold.as_integer_ratio()
    size = len(matrix)
    work = [[b * x - (a if i == j else 0) for j, x in enumerate(row)]
            for i, row in enumerate(matrix)]
    previous = 1
    negative = 0
    for k in range(size):
        row = work[k]
        pivot = row[k]  # Determinant of the leading (k+1)-square block.
        assert pivot != 0
        negative += (pivot > 0) != (previous > 0)
        # Symmetric Bareiss update. The signed LDL pivot is pivot/previous.
        for i in range(k + 1, size):
            target = work[i]
            factor = row[i]
            for j in range(i, size):
                target[j], rem = divmod(target[j] * pivot - factor * row[j],
                                        previous)
                assert rem == 0
        previous = pivot
    return negative


def main():
    kh = knot_complex()
    # At h=-2 there is no spectrum below HIGH. Consequently the low positive
    # eigenvalue at h=-1 must come from d_-1 and also occurs at h=0.
    for h, expected_counts in [(-2, (0, 0)), (-1, (0, 1)), (0, (1, 2))]:
        incoming = kh.differential(h - 1, -1)
        outgoing = kh.differential(h, -1)
        kh.check_d_squared(h - 1, -1)
        betti = incoming.nrows - incoming.rank() - outgoing.rank()
        matrix = laplacian(incoming, outgoing)
        counts = tuple(count_below(matrix, t) for t in (LOW, HIGH))
        assert counts == expected_counts and betti == counts[0]
        print(f"(h,q)=({h},-1) dimension={len(matrix)} Betti={betti} "
              f"counts below the endpoints={counts}")


if __name__ == "__main__":
    main()