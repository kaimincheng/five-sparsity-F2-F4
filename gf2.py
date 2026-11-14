"""Exact binary polynomial arithmetic with nonnegative integer bitmasks.

Bit i is the coefficient of X^i.  Only the Python standard library is used.
"""
from math import gcd


def degree(f):
    return f.bit_length() - 1


def multiply(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        a <<= 1
        b >>= 1
    return out


def divide(a, b):
    """Return quotient, remainder in F_2[X]."""
    assert b > 0
    out = 0
    db = degree(b)
    while a and degree(a) >= db:
        shift = degree(a) - db
        out ^= 1 << shift
        a ^= b << shift
    return out, a


def rem(a, b):
    return divide(a, b)[1]


def gcdpoly(a, b):
    while b:
        a, b = b, rem(a, b)
    return a


def powmod(base, exp, modulus):
    out = 1
    while exp:
        if exp & 1:
            out = rem(multiply(out, base), modulus)
        exp >>= 1
        if exp:
            base = rem(multiply(base, base), modulus)
    return out


def prime_divisors(n):
    result = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            result.append(p)
            while n % p == 0:
                n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        result.append(n)
    return result


def irreducible(f):
    """Rabin's criterion over F_2; degree-one case included."""
    n = degree(f)
    if n < 1 or not (f & 1) or not (f >> n == 1):
        return False
    x = 2
    if powmod(x, 1 << n, f) != rem(x, f):
        return False
    for ell in prime_divisors(n):
        if gcdpoly(f, powmod(x, 1 << (n // ell), f) ^ x) != 1:
            return False
    return True


def order_mod(q, n):
    if n == 1:
        return 1
    assert gcd(q, n) == 1
    r = q % n
    order = 1
    while r != 1:
        r = r * q % n
        order += 1
    return order


def from_exponents(exponents):
    out = 0
    for exponent in exponents:
        out ^= 1 << exponent
    return out


def factor_squarefree_nth_root_unity(n):
    """Factor X^n+1 (n odd) by exhaustive irreducible search, n<=73.

    Degrees cannot exceed ord_n(2), so this requires at most 512
    degree-nine candidates in our applications. All factors occur once.
    """
    assert n > 0 and n % 2 == 1 and n <= 73
    bound = order_mod(2, n)
    factors = []
    target = (1 << n) | 1
    for d in range(1, bound + 1):
        for middle in range(1 << (d - 1)):
            f = (1 << d) | (middle << 1) | 1
            if rem(target, f) == 0 and irreducible(f):
                factors.append(f)
    product = 1
    for f in factors:
        product = multiply(product, f)
    assert product == target, (n, len(factors), bound)
    return factors
