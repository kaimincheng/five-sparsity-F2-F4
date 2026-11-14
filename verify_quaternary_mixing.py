"""Finite check (C8): all sixteen F4 mixing certificates over F2.

Only Python standard-library integer bit arithmetic is used.
"""
from functools import lru_cache
from gf2 import divide, rem

CASES = [
    (35, (1,2,3,4), 0), (65, (1,2,3,4), 0),
    (155, (1,2,3,4), 0), (365, (1,2,3,4), 3),
    (91, (1,2,3,4), 4), (119, (1,2,3,4), 0),
    (217, (1,2,5,6), 0), (511, (1,2,3,4), 0),
    (221, (1,2,3,4), 1), (403, (1,2,3,5), 5),
    (949, (1,2,3,4), 0), (527, (1,2,3,4), 0),
    (1241, (1,2,3,4), 0), (2263, (1,2,3,5), 8),
    (117, (1,2,3,5), 8), (657, (1,2,3,4), 1),
]


def proper_divisors(n):
    return [d for d in range(1, n) if n % d == 0]


@lru_cache(None)
def phi_mod2(n):
    """Cyclotomic Phi_n(X) reduced mod2 by exact polynomial division."""
    work = (1 << n) | 1  # X^n-1 = X^n+1 in characteristic 2
    for d in proper_divisors(n):
        q, r = divide(work, phi_mod2(d))
        assert r == 0, (n, d)
        work = q
    return work


def orbit4(n):
    orbit = [1]
    while (value := orbit[-1] * 4 % n) != 1:
        assert value not in orbit
        orbit.append(value)
    return orbit


def cyc_shift(bits, k, n, mask):
    k %= n
    if k == 0:
        return bits
    return ((bits << k) & mask) | (bits >> (n-k))


def cyc_mul(a,b,n,mask):
    if a.bit_count() > b.bit_count():
        a,b=b,a
    out=0
    while a:
        low=a & -a
        j=low.bit_length()-1
        out ^= cyc_shift(b,j,n,mask)
        a ^= low
    return out


def verify():
    assert len(CASES) == 16
    for d, J, k in CASES:
        orbit = orbit4(d)
        N = max(J)
        mask = (1 << d) - 1
        E = [0] * (N+1)
        E[0] = 1
        for u in orbit:
            for j in range(N,0,-1):
                E[j] ^= cyc_shift(E[j-1],u,d,mask)
        product = 1
        for j in J:
            product = cyc_mul(product,E[j],d,mask)
        remainder = rem(product,phi_mod2(d))
        assert ((remainder >> k) & 1) == 1, (d,J,k)
        print(f'C8: d={d}, J={J}, k={k}, deg(Phi_d mod 2)={phi_mod2(d).bit_length()-1}: PASS')


if __name__ == '__main__':
    verify()
