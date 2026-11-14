"""Finite check (C4): 35 exact affine fifth-norm supports.

Two independent methods are compared: a 5x5 permanent in characteristic
2, and a conjugate product over F16. No third-party libraries needed.
The full verified list of 338 (family, r, slope, intercept) rows is
written using --write-data, or validated against data/ otherwise.
"""
import csv
from collections import defaultdict
from itertools import combinations, permutations, product
from pathlib import Path
import sys

# Each pair (a,b) represents a*k+b.
FAMILIES = {
    'I':   ((6,3), (6,2), (0,2)),
    'II':  ((6,3), (6,2), (6,1)),
    'III': ((6,3), (3,2), (3,1)),
    'IV':  ((18,9), (12,6), (6,3)),
    'V':   ((6,5), (0,4), (0,3)),
    'VI':  ((12,5), (9,4), (3,2)),
    'VII': ((18,1), (12,1), (6,1)),
}
EXPECTED = {
    'I':[9,8,9,9,9], 'II':[9,8,9,9,9],
    'III':[9,9,9,13,8], 'IV':[9,8,17,9,9],
    'V':[11,9,8,9,11], 'VI':[11,9,8,13,9],
    'VII':[17,8,9,9,9],
}
PATH = Path(__file__).resolve().parent / 'data' / 'norm5_affine_supports.csv'


def gf16mul(a, b):
    x = 0
    for i in range(4):
        for j in range(4):
            if (a >> i & 1) and (b >> j & 1):
                x ^= 1 << (i+j)
    for i in range(6, 3, -1):
        if x >> i & 1:
            x ^= 0b10011 << (i-4)
    return x


def gf16pow(a, power):
    out = 1
    for _ in range(power):
        out = gf16mul(out, a)
    return out


ZETA = [gf16pow(2, 3*j) for j in range(5)]
assert ZETA[0] == 1 and ZETA[1] != 1 and gf16mul(ZETA[1],ZETA[4]) == 1


def input_exponents(family, residue):
    # Each exponent has the shape slope*K+intercept, k=5K+residue.
    return [(0,0), (0,1)] + [
        (5*a, a*residue+b) for a,b in FAMILIES[family][::-1]
    ]


def by_conjugates(family, residue):
    exps = input_exponents(family, residue)
    terms = defaultdict(int)
    for indices in product(range(5), repeat=5):
        s = sum(exps[ix][0] for ix in indices)
        c = sum(exps[ix][1] for ix in indices)
        z = sum(j*exps[ix][1] for j,ix in enumerate(indices)) % 5
        terms[(s,c)] ^= ZETA[z]
    result = set()
    for (s,c), coeff in terms.items():
        if coeff:
            assert s%5 == 0 and c%5 == 0 and coeff == 1
            result.add((s//5,c//5))
    return result


def by_determinant(family, residue):
    exps = input_exponents(family, residue)
    matrix = [[set() for j in range(5)] for ell in range(5)]
    for j in range(5):
        for a,b in exps:
            ell = (b+j)%5
            assert a%5 == 0
            q = (a//5, (b+j-ell)//5)
            # Same entry may contain duplicate monomials; work mod 2.
            if q in matrix[ell][j]:
                matrix[ell][j].remove(q)
            else:
                matrix[ell][j].add(q)
    terms = set()
    for rows in permutations(range(5)):
        factors = [matrix[rows[j]][j] for j in range(5)]
        if any(not f for f in factors):
            continue
        for tuple_terms in product(*factors):
            pair = (sum(p[0] for p in tuple_terms),
                    sum(p[1] for p in tuple_terms))
            if pair in terms:
                terms.remove(pair)
            else:
                terms.add(pair)
    return terms


def verify():
    rows = []
    for family in FAMILIES:
        weights = []
        for r in range(5):
            a = by_conjugates(family, r)
            b = by_determinant(family, r)
            assert a == b, (family,r,len(a),len(b),a^b)
            for (s,c),(s1,c1) in combinations(sorted(a),2):
                if s == s1:
                    assert c != c1
                elif (c1-c)%(s-s1)==0:
                    K = (c1-c)//(s-s1)
                    assert K<1, (family,r,K,(s,c),(s1,c1))
            weights.append(len(a))
            for slope,intercept in sorted(a):
                rows.append((family,r,slope,intercept))
        assert weights == EXPECTED[family], (family,weights)
        print(f'C4: {family}, r=0,...,4: weights {weights}; determinant=conjugates: PASS')
    assert len(rows) == 338
    header = ['family','r','slope','intercept']
    if '--write-data' in sys.argv:
        PATH.parent.mkdir(parents=True, exist_ok=True)
        with PATH.open('w',newline='',encoding='utf-8') as fp:
            w = csv.writer(fp)
            w.writerow(header)
            w.writerows(rows)
        print(f'C4: wrote verified data to {PATH}')
    else:
        assert PATH.exists(), f'Missing data file {PATH}; run with --write-data'
        with PATH.open(newline='',encoding='utf-8') as fp:
            read = csv.reader(fp)
            assert next(read) == header
            data = [(x[0],)+tuple(map(int,x[1:])) for x in read]
        assert data == rows, 'CSV file does not match computed supports'
        print('C4: verified all 338 CSV exponent pairs against both algorithms: PASS')


if __name__ == '__main__':
    verify()
