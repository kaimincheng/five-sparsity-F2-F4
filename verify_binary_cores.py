"""Finite check (C9): factor small odd-order cores over F2.

Uses exhaustive enumeration of monic irreducible binary polynomials of
relevant degrees <=9, followed by an exact product identity.
"""
from collections import Counter
from gf2 import factor_squarefree_nth_root_unity, degree

CORES = [3, 5, 7, 9, 15, 21, 31, 63, 73]
EXPECTED_MAX = {3:3, 5:5, 7:3, 9:3, 15:5, 21:5, 31:5, 63:5, 73:5}
EXPECTED_NUMBER = {3:2, 5:2, 7:3, 9:3, 15:5, 21:6, 31:7, 63:13, 73:9}


def verify():
    for n in CORES:
        fac = factor_squarefree_nth_root_unity(n)
        weights = Counter(f.bit_count() for f in fac)
        assert len(fac) == EXPECTED_NUMBER[n], (n, len(fac))
        assert max(weights) == EXPECTED_MAX[n], (n, weights)
        assert all(f.bit_count() <= 5 for f in fac)
        if n == 73:
            d9 = [f for f in fac if degree(f) == 9]
            assert len(d9) == 8
            assert all(f.bit_count() in (3, 5) for f in d9)
        if n == 21:
            f1 = sum(1 << e for e in [6,4,2,1,0])
            f2 = sum(1 << e for e in [6,5,4,2,0])
            assert f1 in fac and f2 in fac
        print(f'C9: n={n}, irreducible factors={len(fac)}, max weight={max(weights)}, weight histogram={dict(sorted(weights.items()))}: PASS')


if __name__ == '__main__':
    verify()
