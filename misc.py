
from math import gcd, sqrt
import itertools
from typing import Iterator

tiny_primes = 2, 3, 5, 7, 11, 13, 17, 19
small_primes = frozenset(itertools.chain(
    tiny_primes,
    (x for x in range(23, 500) if not any(x % p == 0 for p in tiny_primes))))

modest_primes = frozenset(itertools.chain(
    small_primes,
    (x for x in range(501, 250000, 2)
     if not any(x % p == 0 for p in small_primes))))

def jacobi(a: int, n: int) -> int:
    assert n > 0 and n & 1 == 1         # Must be positive, odd.
    sign = 1
    if n == 1:
        return 1
    while True:
        a %= n
        if a == 0:
            return 0
        # Divide out twos.  (2|n) = 1 if n ≡ ±1 mod 8 and -1 if n ≡ ±3 mod 8.
        while a & 3 == 0:
            a >>= 2
        if a & 1 == 0:
            a >>= 1
            if n & 7 in (3,5):
                sign = -sign

        if a == 1:
            return sign

        # Apply reciprocity.
        if 3 & a & n == 3:
            sign = -sign

        a, n = n, a

def is_square(a: int) -> bool:
    if a < 0:
        return False
    if a < 2:
        return True
    # Find the smallest power of two such that s*s >= a.
    s = 2
    ssq = 4
    lower = 1
    while ssq < a:
        lower = s
        s = 2 * s
        ssq = 4 * ssq
    # lower² < a
    # Hence if we find a < s² < a + 2lower + 1 then we know that s is the only
    # square in that range.
    stop = a + 2 * lower + 1
    while ssq >= stop:
        # This gives s² ≥ a if the division is not floored, i.e., s ≥ √a before
        # flooring and so ⌊s⌋ ≥ ⌊√a⌋.
        s = (a + ssq) // (2 * s)
        ssq = s * s
    return ssq == a

if __name__ == '__main__':
    assert list(tiny_primes) == sorted(tiny_primes)
    assert list(small_primes) == sorted(small_primes)

    for i in range(2000*0):
        assert is_square(i*i)
        for j in range(i*i + 1, (i+1)*(i+1)):
            assert not is_square(j)

    import cmath
    def rc(c):
        s = 1048576
        return round(c.real*s)/s + round(c.imag*s)/s*1j
    for p in range(3,500,2):
        if p == 2:
            continue
        S = 0j
        for k in range(1, p):
            S += jacobi(k, p) * cmath.rect(1, 2 * k * cmath.pi / p)
        print(p, rc(S*S), rc(S / sqrt(p)))
