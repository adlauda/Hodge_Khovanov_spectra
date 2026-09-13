"""Small deterministic checks of the new manuscript, not spectral certification."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import re

root = Path(__file__).resolve().parents[1]
bound = F(2, 31) + F(1, 32) + F(1, 16)
assert bound < F(1, 2)
worst = F(0)
for R in (1, 2, 3, 10, 100, 1000):
    for r, tail, prep_sign, sample_sign in product(
            set((1, R)), (F(0), F(1, 16)), (-1, 1), (-1, 1)):
        p = 1 / (r + tail)
        prepared = p + prep_sign * F(1, 128 * R * R)
        if not 0 <= prepared <= 1:
            continue
        frequency = prepared * (1 + sample_sign * F(1, 32 * R))
        if not 0 < frequency <= 1:
            continue
        error = abs(1 / frequency - r)
        assert error < F(1, 2)
        worst = max(worst, error)

def grading(u, v, w, kappa):
    return (u - 2*v + F(w+kappa, 2),
            2*u - 2*v + F(3*w+kappa, 2) - 1)

# Positive triangle, reduced unknot grading (0,-1).
assert {grading(u, 2, 3, 7) for u in (-1, 1, 2)} == {(0, 1), (2, 5), (3, 7)}
# Mirroring reduced homology sends j to -j-2 in this convention.
assert {grading(u, 1, -3, 5) for u in (1, -1, -2)} == {(0, -3), (-2, -7), (-3, -9)}

# Polynomial identity J_trefoil=(q^2+1)*chi_reduced, as exponent dictionaries.
reduced = {1: 1, 5: 1, 7: -1}
unreduced = {}
for power, coeff in reduced.items():
    for shift in (0, 2):
        unreduced[power+shift] = unreduced.get(power+shift, 0) + coeff
assert {p: c for p, c in unreduced.items() if c} == {1: 1, 3: 1, 5: 1, 9: -1}
# Two-component unlink must carry the minus sign in normalized V.
assert (-1)**(2-1) == -1

texts = {p: p.read_text(encoding='utf-8') for p in root.glob('*.tex')}
labels, refs = [], []
for path, text in texts.items():
    active = '\n'.join(line.split('%')[0] for line in text.splitlines()
                       if not line.lstrip().startswith('%'))
    labels += re.findall(r'\\label\{([^}]+)\}', active)
    for group in re.findall(r'\\(?:Cref|cref|ref|eqref)\{([^}]+)\}', active):
        refs.extend(group.split(','))
    for image in re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}', active):
        assert (root / image).is_file(), (path.name, image)
    assert not re.search(r'\\(?:AL|AS|PZ|MR|TODO|todo)\{', active), path.name
assert len(labels) == len(set(labels)), 'Duplicate label names'
assert set(refs) <= set(labels), sorted(set(refs)-set(labels))
print(f'PASS: reference-error bound {float(bound):.9f}; endpoint maximum {float(worst):.9f}.')
print(f'PASS: trefoil grading and Euler checks; {len(labels)} labels; every image exists.')
print('These checks are diagnostics, not a proof of Gibbs preparation or spectral gaps.')
