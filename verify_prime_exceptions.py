"""Finite checks (C5) and (C6) of the paper. Standard library only."""
from gf2 import from_exponents, irreducible, order_mod, powmod, prime_divisors


def verify_small_orders():
    expected = {3: [7], 5: [31], 7: [127], 9: [73],
                11: [23, 89], 13: [8191]}
    for t, primes in expected.items():
        n = (1 << t) - 1
        actual = [p for p in prime_divisors(n) if order_mod(2, p) == t]
        assert actual == sorted(primes), (t, actual, primes)
        print(f'C5: ord_p(2)={t}: primes {actual}')


def verify_binary_factors():
    data = {
        127: [0, 2, 3, 4, 5, 6, 7],
        23: [0, 1, 5, 6, 7, 9, 11],
        89: [0, 1, 2, 3, 5, 6, 7, 8, 11],
    }
    for p, es in data.items():
        f = from_exponents(es)
        assert irreducible(f)
        assert f.bit_length() - 1 == order_mod(2, p)
        assert powmod(2, p, f) == 1
        assert all(powmod(2, p // z, f) != 1 for z in prime_divisors(p))
        assert len(es) > 5
        print(f'C5: p={p}, degree={max(es)}, weight={len(es)}: PASS')


def gf4_mul(a, b):
    """Bits encode F_4=F_2[w]/(w^2+w+1), w=2, w^2=3."""
    out = 0
    for i in range(2):
        for j in range(2):
            if (a >> i) & 1 and (b >> j) & 1:
                out ^= 1 << (i + j)
    if out & 4:
        out ^= 7
    return out


def gf4_poly_mul(a, b):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            e = i + j
            out[e] = out.get(e, 0) ^ gf4_mul(x, y)
    return {e: c for e, c in out.items() if c}


def verify_phi13():
    f = {6: 1, 5: 2, 3: 3, 1: 2, 0: 1}
    g = {6: 1, 5: 3, 3: 2, 1: 3, 0: 1}
    product = gf4_poly_mul(f, g)
    assert product == {i: 1 for i in range(13)}
    assert order_mod(4, 13) == 6
    print('C6: Phi_13 factorization into two degree-six F4 factors: PASS')


if __name__ == '__main__':
    verify_small_orders()
    verify_binary_factors()
    verify_phi13()
