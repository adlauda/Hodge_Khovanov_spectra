"""Independent, exact rational Khovanov cube engine for diagram stress tests.

This module does not import or use the manuscript's verifier.  It constructs
every required resolution directly from the external edges of a closed diagram,
uses Q[X]/(X^2), and performs sparse Gaussian elimination with fractions.Fraction.

Interface
---------
    d = Diagram(crossings, edges, npos=None, nneg=0, marks={"p": vertex})
    kh = Khovanov(d)
    kh.chain(i, j)              # ChainGroup: .basis, .index, len(...)
    kh.differential(i, j)       # SparseMatrix: .columns, .nrows, .ncols, .rank()
    kh.homology(i, j)           # Homology: .dimension, .representatives, .project(v)
    kh.dot_map(i, j, "p")       # H^{i,j} -> H^{i,j-2}, exact sparse matrix
    kh.check_d_squared(i, j)    # raises if the two differentials do not compose to 0
    kh.check_dot_chain_map(i,j,"p")  # raises if d X_p != X_p d in that bidegree

Each crossing is (NW, NE, SE, SW), with distinct hashable vertex names.  Its
0-smoothing pairs (NW, SW), (NE, SE); its 1-smoothing pairs (NW, NE), (SE, SW).
External edges join half-edges outside crossing discs.  A crossing vertex has
exactly one incident external edge, and each subdivision vertex has exactly two.
A crossingless circle may be represented by a loop edge (v,v).  This data must
describe an ordinary link diagram for the usual Khovanov interpretation; a
nonplanar abstract matching with a 1-to-1 saddle is rejected.

All crossings use the indicated local smoothing convention.  The npos/nneg
counts change normalization only; they do not change local smoothings.  Thus
geometrically negative crossings must be encoded with the appropriate rotated
local half-edge order by the caller.  For a state with r 1-resolutions and k
circles, a label bit equal to 1 denotes X (degree -1), and 0 denotes 1 (degree +1):
    i = r - nneg
    j = k - 2 * number_of_X_labels + r + npos - 2*nneg.
The cube-edge sign is (-1)^(number of 1-bits preceding the changed crossing).

Only the states and chain groups requested are generated.  Top-row or narrow
quantum-degree calculations at 12--16 crossings are generally modest; computing
all middle chain groups can be exponentially expensive.  No modular arithmetic,
floating-point arithmetic, or unproved rank shortcuts are used.

Run this file directly for independent regression tests: the unknot, unlink,
positive Hopf link, and right-handed trefoil, including d^2 and marked dot maps.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from typing import Dict, Hashable, Iterable, List, Mapping, Optional, Sequence, Tuple


Vector = Dict[int, Fraction]


def popcount(n: int) -> int:
    return bin(n).count("1")


def _add_scaled(target: Vector, source: Mapping[int, Fraction], scale) -> None:
    """target += scale * source, removing zero entries."""
    if not scale:
        return
    for row, value in source.items():
        new = target.get(row, 0) + scale * value
        if new:
            target[row] = new
        else:
            target.pop(row, None)


def _fraction_vector(v: Mapping[int, object]) -> Vector:
    return {k: Fraction(x) for k, x in v.items() if x}


class Echelon:
    """Sparse column echelon basis, with exact optional companion vectors."""

    def __init__(self):
        self.pivots: Dict[int, Vector] = {}
        self.companions: Dict[int, Vector] = {}
        self.tags: Dict[int, object] = {}

    def reduce(self, vector, companion=None, coordinates=None):
        v = _fraction_vector(vector)
        c = _fraction_vector(companion) if companion is not None else None
        while v:
            pivot = min(v)
            if pivot not in self.pivots:
                break
            factor = v[pivot]
            _add_scaled(v, self.pivots[pivot], -factor)
            if c is not None:
                _add_scaled(c, self.companions[pivot], -factor)
            if coordinates is not None:
                tag = self.tags[pivot]
                if tag is not None:
                    coordinates[tag] = coordinates.get(tag, 0) + factor
                    if not coordinates[tag]:
                        del coordinates[tag]
        return v, c

    def insert(self, vector, companion=None, tag=None):
        """Return (independent, reduced_vector, reduced_companion)."""
        v, c = self.reduce(vector, companion)
        if not v:
            return False, v, c
        pivot = min(v)
        scale = v[pivot]
        v = {r: x / scale for r, x in v.items()}
        if c is not None:
            c = {r: x / scale for r, x in c.items()}
            self.companions[pivot] = c
        self.pivots[pivot] = v
        self.tags[pivot] = tag
        return True, v, c

    def __len__(self):
        return len(self.pivots)


@dataclass
class SparseMatrix:
    nrows: int
    columns: List[Vector]
    _rank: Optional[int] = field(default=None, repr=False)

    @property
    def ncols(self):
        return len(self.columns)

    def apply(self, vector: Mapping[int, Fraction]) -> Vector:
        result: Vector = {}
        for col, coeff in vector.items():
            _add_scaled(result, self.columns[col], coeff)
        return result

    def rank(self) -> int:
        if self._rank is None:
            echelon = Echelon()
            for col in self.columns:
                echelon.insert(col)
            self._rank = len(echelon)
        return self._rank

    def kernel(self) -> List[Vector]:
        echelon = Echelon()
        result = []
        for col, vector in enumerate(self.columns):
            independent, _, companion = echelon.insert(vector, {col: Fraction(1)})
            if not independent:
                result.append(companion)
        self._rank = len(echelon)
        return result

    def dense(self) -> List[List[Fraction]]:
        result = [[Fraction(0) for _ in self.columns] for _ in range(self.nrows)]
        for col, vector in enumerate(self.columns):
            for row, value in vector.items():
                result[row][col] = value
        return result

    def compose(self, other: "SparseMatrix") -> "SparseMatrix":
        """Return self after other."""
        if self.ncols != other.nrows:
            raise ValueError("Matrix dimensions do not compose")
        return SparseMatrix(self.nrows, [self.apply(col) for col in other.columns])


@dataclass
class Diagram:
    crossings: Sequence[Tuple[Hashable, Hashable, Hashable, Hashable]]
    edges: Sequence[Tuple[Hashable, Hashable]]
    npos: Optional[int] = None
    nneg: int = 0
    marks: Mapping[str, Hashable] = field(default_factory=dict)

    def __post_init__(self):
        self.crossings = tuple(tuple(c) for c in self.crossings)
        self.edges = tuple(tuple(e) for e in self.edges)
        self.marks = dict(self.marks)
        c = len(self.crossings)
        if self.npos is None:
            self.npos = c - self.nneg
        if self.npos < 0 or self.nneg < 0 or self.npos + self.nneg != c:
            raise ValueError("npos and nneg must be nonnegative and sum to crossing count")
        crossing_vertices = []
        for cross in self.crossings:
            if len(cross) != 4:
                raise ValueError("Each crossing requires four half-edge vertices")
            crossing_vertices.extend(cross)
        if len(set(crossing_vertices)) != len(crossing_vertices):
            raise ValueError("Crossing half-edge vertices must all be distinct")
        degrees = {}
        for edge in self.edges:
            if len(edge) != 2:
                raise ValueError("Each external edge requires two endpoints")
            for v in edge:
                degrees[v] = degrees.get(v, 0) + 1
        crossing_set = set(crossing_vertices)
        for v in crossing_set | set(degrees):
            expected = 1 if v in crossing_set else 2
            if degrees.get(v, 0) != expected:
                raise ValueError("External degree of %r is %s; expected %s" %
                                 (v, degrees.get(v, 0), expected))
        if any(v not in degrees for v in self.marks.values()):
            raise ValueError("A marked vertex is absent from the diagram")


@dataclass(frozen=True)
class Resolution:
    circles: Tuple[Tuple[int, ...], ...]
    vertex_circle: Tuple[int, ...]

    @property
    def count(self):
        return len(self.circles)


@dataclass
class ChainGroup:
    i: int
    j: int
    basis: List[Tuple[int, int]]
    index: Dict[Tuple[int, int], int]

    def __len__(self):
        return len(self.basis)


@dataclass
class Homology:
    i: int
    j: int
    chain: ChainGroup
    representatives: List[Vector]
    boundary_rank: int
    outgoing_rank: int
    echelon: Echelon = field(repr=False)

    @property
    def dimension(self):
        return len(self.representatives)

    def project(self, cycle: Mapping[int, Fraction]) -> Vector:
        """Coordinates of a cycle in this homology basis; rejects noncycles."""
        result: Vector = {}
        remainder, _ = self.echelon.reduce(cycle, coordinates=result)
        if remainder:
            raise ValueError("The supplied vector is not in the cycle space")
        return result

    def coordinates(self, cycle) -> Vector:
        """Like project; additionally accept (state,labelmask) keys directly."""
        converted = {}
        for key, value in cycle.items():
            row = self.chain.index[key] if isinstance(key, tuple) else key
            converted[row] = converted.get(row, 0) + value
        return self.project(converted)


class Khovanov:
    def __init__(self, diagram: Diagram):
        self.diagram = diagram
        self.c = len(diagram.crossings)
        # Input appearance order avoids requiring comparable vertex identifiers.
        names = dict.fromkeys(v for edge in diagram.edges for v in edge)
        self.vertices = tuple(names)
        self.vertex_index = {v: n for n, v in enumerate(self.vertices)}
        self.crossings = tuple(tuple(self.vertex_index[v] for v in cross)
                               for cross in diagram.crossings)
        self.edges = tuple((self.vertex_index[a], self.vertex_index[b])
                           for a, b in diagram.edges)
        self.marks = {name: self.vertex_index[v] for name, v in diagram.marks.items()}
        # Caches belong to this engine rather than global decorators.  An engine
        # may therefore be garbage-collected after an enumeration test case.
        for method in ("resolution", "chain", "saddle", "differential", "homology"):
            setattr(self, method, lru_cache(maxsize=None)(getattr(self, method)))

    def resolution(self, state: int) -> Resolution:
        if state < 0 or state >= 1 << self.c:
            raise ValueError("State is outside the cube")
        parent = list(range(len(self.vertices)))
        size = [1] * len(parent)

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        def union(a, b):
            a, b = find(a), find(b)
            if a != b:
                if size[a] < size[b]:
                    a, b = b, a
                parent[b] = a
                size[a] += size[b]

        for a, b in self.edges:
            union(a, b)
        for q, (a, b, c, d) in enumerate(self.crossings):
            if state & (1 << q):
                union(a, b)
                union(c, d)
            else:
                union(a, d)
                union(b, c)
        components = {}
        for v in range(len(parent)):
            components.setdefault(find(v), []).append(v)
        circles = tuple(sorted((tuple(vs) for vs in components.values()), key=lambda vs: vs[0]))
        assignment = [0] * len(parent)
        for circle, vs in enumerate(circles):
            for v in vs:
                assignment[v] = circle
        return Resolution(circles, tuple(assignment))

    def chain(self, i: int, j: int) -> ChainGroup:
        r = i + self.diagram.nneg
        basis = []
        if 0 <= r <= self.c:
            shift = r + self.diagram.npos - 2 * self.diagram.nneg
            for ones in combinations(range(self.c), r):
                state = sum(1 << q for q in ones)
                k = self.resolution(state).count
                twice_x = k + shift - j
                if twice_x % 2 or twice_x < 0 or twice_x > 2 * k:
                    continue
                for x_circles in combinations(range(k), twice_x // 2):
                    labels = sum(1 << q for q in x_circles)
                    basis.append((state, labels))
        return ChainGroup(i, j, basis, {b: n for n, b in enumerate(basis)})

    def saddle(self, state: int, crossing: int):
        """Topology of one 0->1 saddle, derived from actual resolution circles."""
        if state & (1 << crossing):
            raise ValueError("Saddle source already has a 1-resolution")
        before = self.resolution(state)
        after = self.resolution(state | (1 << crossing))
        incidence = [set(after.vertex_circle[v] for v in vs) for vs in before.circles]
        reverse = [set() for _ in after.circles]
        for old, targets in enumerate(incidence):
            for new in targets:
                reverse[new].add(old)
        if after.count == before.count - 1:
            merged = [new for new, sources in enumerate(reverse) if len(sources) == 2]
            if len(merged) != 1 or any(len(targets) != 1 for targets in incidence):
                raise ValueError("Malformed merge saddle")
            new = merged[0]
            old1, old2 = sorted(reverse[new])
            common = tuple((old, next(iter(targets))) for old, targets in enumerate(incidence)
                           if old != old1 and old != old2)
            return "merge", (old1, old2), (new,), common
        if after.count == before.count + 1:
            split = [old for old, targets in enumerate(incidence) if len(targets) == 2]
            if len(split) != 1 or any(len(sources) != 1 for sources in reverse):
                raise ValueError("Malformed split saddle")
            old = split[0]
            new1, new2 = sorted(incidence[old])
            common = tuple((other, next(iter(targets))) for other, targets in enumerate(incidence)
                           if other != old)
            return "split", (old,), (new1, new2), common
        raise ValueError("A saddle does not change circle count by one; check diagram planarity")

    def _saddle_labels(self, state, labels, crossing):
        kind, old, new, common = self.saddle(state, crossing)
        result = 0
        for a, b in common:
            if labels & (1 << a):
                result |= 1 << b
        if kind == "merge":
            x_count = bool(labels & (1 << old[0])) + bool(labels & (1 << old[1]))
            if x_count == 2:
                return ()
            if x_count == 1:
                result |= 1 << new[0]
            return (result,)
        if labels & (1 << old[0]):
            return (result | (1 << new[0]) | (1 << new[1]),)
        return (result | (1 << new[0]), result | (1 << new[1]))

    def differential(self, i: int, j: int) -> SparseMatrix:
        source = self.chain(i, j)
        target = self.chain(i + 1, j)
        columns = []
        for state, labels in source.basis:
            vector = {}
            for q in range(self.c):
                if state & (1 << q):
                    continue
                sign = -1 if popcount(state & ((1 << q) - 1)) % 2 else 1
                target_state = state | (1 << q)
                for target_labels in self._saddle_labels(state, labels, q):
                    row = target.index[(target_state, target_labels)]
                    vector[row] = vector.get(row, 0) + sign
            columns.append({r: Fraction(x) for r, x in vector.items() if x})
        return SparseMatrix(len(target), columns)

    def homology(self, i: int, j: int) -> Homology:
        chain = self.chain(i, j)
        incoming = self.differential(i - 1, j)
        outgoing = self.differential(i, j)
        quotient = Echelon()
        for boundary in incoming.columns:
            if outgoing.apply(boundary):
                raise AssertionError("d^2 is not zero at (%s,%s)" % (i - 1, j))
            quotient.insert(boundary)
        boundary_rank = len(quotient)
        incoming._rank = boundary_rank
        representatives = []
        for cycle in outgoing.kernel():
            independent, reduced, _ = quotient.insert(cycle, tag=len(representatives))
            if independent:
                representatives.append(reduced)
        result = Homology(i, j, chain, representatives, boundary_rank, outgoing.rank(), quotient)
        expected = len(chain) - boundary_rank - outgoing.rank()
        if result.dimension != expected:
            raise AssertionError("Homology dimension and rank-nullity disagree")
        return result

    def dot_chain(self, i: int, j: int, mark: str) -> SparseMatrix:
        vertex = self.marks[mark]
        source, target = self.chain(i, j), self.chain(i, j - 2)
        columns = []
        for state, labels in source.basis:
            circle = self.resolution(state).vertex_circle[vertex]
            if labels & (1 << circle):
                columns.append({})
            else:
                columns.append({target.index[(state, labels | (1 << circle))]: Fraction(1)})
        return SparseMatrix(len(target), columns)

    def dot_map(self, i: int, j: int, mark: str) -> SparseMatrix:
        source, target = self.homology(i, j), self.homology(i, j - 2)
        chain_map = self.dot_chain(i, j, mark)
        return SparseMatrix(target.dimension,
                            [target.project(chain_map.apply(rep)) for rep in source.representatives])

    def check_d_squared(self, i: int, j: int) -> bool:
        a, b = self.differential(i, j), self.differential(i + 1, j)
        for col, image in enumerate(a.columns):
            if b.apply(image):
                raise AssertionError("d^2 != 0 at (%s,%s), column %s" % (i, j, col))
        return True

    def check_dot_chain_map(self, i: int, j: int, mark: str) -> bool:
        d = self.differential(i, j)
        dot_after = self.dot_chain(i + 1, j, mark)
        dot_before = self.dot_chain(i, j, mark)
        d_after = self.differential(i, j - 2)
        for col in range(d.ncols):
            if dot_after.apply(d.columns[col]) != d_after.apply(dot_before.columns[col]):
                raise AssertionError("Dot fails chain-map identity at (%s,%s), column %s" % (i,j,col))
        return True

    def support(self) -> Dict[Tuple[int, int], int]:
        """Full bigraded rational ranks; intended only for small diagrams."""
        degrees = set()
        for state in range(1 << self.c):
            r = popcount(state)
            k = self.resolution(state).count
            shift = r + self.diagram.npos - 2 * self.diagram.nneg
            for x_count in range(k + 1):
                degrees.add((r - self.diagram.nneg, k - 2*x_count + shift))
        return {(i,j): self.homology(i,j).dimension for i,j in sorted(degrees)
                if self.homology(i,j).dimension}


def two_strand_braid_closure(crossing_count: int) -> Diagram:
    """Standard closure of the positive two-strand braid sigma_1^m."""
    if crossing_count == 0:
        return Diagram([], [(0,0), (1,1)], marks={"left":0,"right":1})
    crossings = [tuple(4*q + p for p in range(4)) for q in range(crossing_count)]
    edges = []
    for q, cross in enumerate(crossings):
        nxt = crossings[(q+1) % crossing_count]
        edges.extend([(cross[3], nxt[0]), (cross[2], nxt[1])])
    return Diagram(crossings, edges, marks={"left":crossings[0][0], "right":crossings[0][1]})


def braid_closure(strands: int, word: Sequence[int]) -> Diagram:
    """Close a signed braid word (one-based generators); helper for validation.

    A negative crossing is encoded by a quarter-turn of its local half-edge
    order, interchanging its two smoothings.  All braid strands are oriented
    downwards, so the normalization is the number of positive/negative letters.
    """
    if strands < 1:
        raise ValueError("A braid must have at least one strand")
    tops = [("top", p) for p in range(strands)]
    active = list(tops)
    edges, crossings = [], []
    for q, generator in enumerate(word):
        p = abs(generator) - 1
        if generator == 0 or p + 1 >= strands:
            raise ValueError("Braid generator is out of range")
        nw, ne, se, sw = [(q, a) for a in range(4)]
        edges.extend([(active[p], nw), (active[p+1], ne)])
        active[p], active[p+1] = sw, se
        crossings.append((nw,ne,se,sw) if generator > 0 else (ne,se,sw,nw))
    edges.extend(zip(active, tops))
    return Diagram(crossings, edges, nneg=sum(g < 0 for g in word),
                   marks={str(p):v for p,v in enumerate(tops)})


def run_self_tests():
    # A separate dense reduced-row-echelon implementation checks the sparse
    # column algorithm, including rectangular matrices and singular kernels.
    import random
    rng = random.Random(173421)

    def dense_rank(rows):
        rows = [[Fraction(x) for x in row] for row in rows]
        m, n, pivot = len(rows), len(rows[0]), 0
        for col in range(n):
            selected = next((r for r in range(pivot,m) if rows[r][col]), None)
            if selected is None:
                continue
            rows[pivot], rows[selected] = rows[selected], rows[pivot]
            factor = rows[pivot][col]
            rows[pivot] = [x/factor for x in rows[pivot]]
            for r in range(m):
                if r != pivot:
                    factor = rows[r][col]
                    rows[r] = [x-factor*y for x,y in zip(rows[r],rows[pivot])]
            pivot += 1
            if pivot == m:
                break
        return pivot

    for m in range(1,9):
        for n in range(1,9):
            for _ in range(5):
                rows = [[rng.randrange(-3,4) for _ in range(n)] for _ in range(m)]
                matrix = SparseMatrix(m, [{r:Fraction(rows[r][c]) for r in range(m) if rows[r][c]}
                                          for c in range(n)])
                rank = dense_rank(rows)
                kernel = matrix.kernel()
                assert matrix.rank() == rank and len(kernel) == n-rank
                assert all(not matrix.apply(v) for v in kernel)
                assert SparseMatrix(n,kernel).rank() == len(kernel)
    print("PASS 320 randomized exact rank/kernel comparisons against separate dense RREF")
    cases = [
        ("unknot", Diagram([], [(0,0)], marks={"p":0}), {(0,-1):1,(0,1):1}),
        ("two-component unlink", two_strand_braid_closure(0), {(0,-2):1,(0,0):2,(0,2):1}),
        ("positive kink", two_strand_braid_closure(1), {(0,-1):1,(0,1):1}),
        ("positive Hopf", two_strand_braid_closure(2), {(0,0):1,(0,2):1,(2,4):1,(2,6):1}),
        ("right-handed trefoil", two_strand_braid_closure(3), {(0,1):1,(0,3):1,(2,5):1,(3,9):1}),
        ("negative kink", braid_closure(2,[-1]), {(0,-1):1,(0,1):1}),
        ("Reidemeister II", braid_closure(2,[1,-1]), {(0,-2):1,(0,0):2,(0,2):1}),
        ("left-handed trefoil", braid_closure(2,[-1,-1,-1]),
         {(-3,-9):1,(-2,-5):1,(0,-3):1,(0,-1):1}),
    ]
    for name, diagram, expected in cases:
        kh = Khovanov(diagram)
        actual = kh.support()
        assert actual == expected, (name, actual, expected)
        for i in range(-diagram.nneg-1, diagram.npos+1):
            for j in range(-3*kh.c-4, 3*kh.c+5):
                kh.check_d_squared(i,j)
                for mark in kh.marks:
                    kh.check_dot_chain_map(i,j,mark)
        for (i,j) in actual:
            for mark in kh.marks:
                dot = kh.dot_map(i,j,mark)
                second = kh.dot_map(i,j-2,mark)
                assert all(not second.apply(col) for col in dot.columns), (name,"X^2",i,j,mark)
        print("PASS", name, actual)
    kh = Khovanov(cases[0][1])
    assert kh.dot_map(0,1,"p").rank() == 1
    kh = Khovanov(two_strand_braid_closure(2))
    # On Hopf top homology the two component dot actions agree; on the lower
    # row they differ by a sign.  This catches circle tracking, basis projection,
    # and the characteristic-zero Frobenius algebra simultaneously.
    top_l, top_r = kh.dot_map(2,6,"left"), kh.dot_map(2,6,"right")
    assert top_l.rank() == 1 and top_l.columns == top_r.columns
    low_l, low_r = kh.dot_map(0,2,"left"), kh.dot_map(0,2,"right")
    assert low_l.rank() == 1
    assert low_l.columns == [{r:-v for r,v in col.items()} for col in low_r.columns]
    print("PASS exact marked dot actions, chain-map identities, and nilpotence")
    left = Khovanov(braid_closure(3,[1,2,1,2,2]))
    right = Khovanov(braid_closure(3,[2,1,2,2,2]))
    assert left.support() == right.support()
    print("PASS Reidemeister III with a nontrivial exterior braid")
    return True


if __name__ == "__main__":
    run_self_tests()
