"""Finite check (C7). Verify all eight F2 mixed-order witnesses."""
from gf2 import (from_exponents, irreducible, order_mod, powmod,
                 prime_divisors)

DATA = {
    35:   [0,1,2,4,7,8,9,10,12],
    93:   [0,1,2,5,7,9,10],
    155:  [0,1,3,6,11,12,14,15,18,19,20],
    217:  [0,2,8,10,13,14,15],
    219:  [0,4,6,8,10,13,14,16,18],
    365:  [0,8,9,16,18,24,26,27,32,33,34,35,36],
    511:  [0,1,2,3,4,8,9],
    2263: [0,9,16,18,26,34,35,36,40,42,45],
}


def verify():
    assert len(DATA) == 8
    for d, exponents in DATA.items():
        f = from_exponents(exponents)
        deg = max(exponents)
        assert deg == order_mod(2, d)
        assert len(exponents) > 5
        assert irreducible(f), d
        assert powmod(2, d, f) == 1, d
        for ell in prime_divisors(d):
            assert powmod(2, d // ell, f) != 1, (d, ell)
        print(f'C7: d={d}, degree={deg}, weight={len(exponents)}: PASS')


if __name__ == '__main__':
    verify()
