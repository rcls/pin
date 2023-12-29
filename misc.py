
from math import gcd, sqrt
import itertools
from typing import Any, Iterator, Tuple

class FoundFactor(Exception): pass

def sieve_primes_to(n: int) -> Iterator[int]:
    yield 2
    n2 = n // 2
    # 2n+1 is at index n
    sieve = bytearray(n2)
    sieve[0] = 1
    for j in range(n2):
        if sieve[j]:
            continue
        p = j * 2 + 1
        yield p
        for k in range(p * p - 1 >> 1, n2, p):
            sieve[k] = 1

tiny_prime_limit = 22
tiny_primes = tuple(sieve_primes_to(tiny_prime_limit))
small_prime_limit = 500
small_primes = frozenset(sieve_primes_to(small_prime_limit))
modest_prime_limit = 250000
modest_primes_list = tuple(sieve_primes_to(modest_prime_limit))
modest_primes = set(modest_primes_list)

def split_twos(n: int) -> Tuple[int, int]:
    d = n
    s = 0
    while d & 3 == 0:
        d >>= 2
        s += 2
    if d & 1 == 0:
        d >>= 1
        s += 1
    return d, s

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
    # Find the smallest power of two such that s² >= a.
    #
    # If a is a power of 2, a = 2^{a.bit_length() - 1} and so take
    # 2^{a.bit_length()/2}.
    #
    if a.bit_length() == 1:
        twos = a.bit_length() // 2
    else:
        twos = (a.bit_length() + 1) // 2
    s = 1 << twos
    ssq = 1 << 2*twos
    assert ssq >= a
    lower = s // 2
    assert lower * lower <= a
    # s/2 ≤ √a, so if we find s₁ with √a² = a ≤ s₁² < a + s + 1 ≤ a + 2√a + 1 =
    # (√a+1)² then we know that s₁² is the only square in that range.
    stop = a + s + 1
    while ssq >= stop:
        # This would give s² ≥ a if the division were not floored, i.e., s ≥ √a
        # before flooring and so ⌊s⌋ ≥ ⌊√a⌋.
        s = (a // s + s) // 2
        ssq = s * s
    return ssq == a

def floor_sqrt(a: int) -> Tuple[int, int]:
    '''Return (⌊√a⌋, ⌊√a⌋²).'''
    if a < 2:
        assert a >= 0
        return a, a
    # Find the smallest power of two such that s² >= a.
    twos = (a.bit_length() + 1) // 2
    s = 1 << twos
    ssq = 1 << 2*twos
    while ssq > a:
        # This gives s² ≥ a if the division is not floored, i.e., s ≥ √a before
        # flooring and so ⌊s⌋ ≥ ⌊√a⌋.  At the last step we will start with the
        # previous value just above √a and the flooring will give ≤ √a.
        s = (a // s + s) // 2
        ssq = s * s
    return s, ssq

def floor_cbrt(a: int) -> Tuple[int, int]:
    '''Return (⌊∛a⌋, ⌊∛a⌋³).'''
    if a < 8:
        assert a >= 0                   # Lazy.
        return (0, 0) if a == 0 else (1, 1)
    c = 1 << ((a.bit_length() + 1) // 3)
    # Do at least one step, in case we start with c < ⌊∛a⌋.
    while True:
        c = (a // (c*c) + 2 * c) // 3
        # (·)³ is increasing with positive second derivative, so the Newton step
        # always over-estimates (before flooring).  Hence we now have c ≥ ⌊∛a⌋.
        # When we get to also c ≤ ∛a, i.e., c³ ≤ a, we are done.
        cube_c = c * c * c
        if cube_c <= a:
            return c, cube_c

def euclid(a: int, b: int) -> Tuple[int, int, int]:
    if abs(a) < abs(b):
        x, y = a, b
        x_a, x_b = 1, 0
        y_a, y_b = 0, 1
    else:
        x, y = b, a
        x_a, x_b = 0, 1
        y_a, y_b = 1, 0
    while x:
        q = y // x
        y = y % x
        y_a = y_a - q * x_a
        y_b = y_b - q * x_b
        x, y = y, x
        x_a, y_a = y_a, x_a
        x_b, y_b = y_b, x_b
        assert x == x_a * a + x_b * b
        assert y == y_a * a + y_b * b
    return y_a, y_b, y

def test_tiny_primes() -> None:
    assert list(tiny_primes) == sorted(tiny_primes)
    assert not 0 in tiny_primes
    assert not 1 in tiny_primes
    for i in range(2, tiny_prime_limit):
        assert (i in tiny_primes) == all(i % f != 0 for f in range(2, i))

def test_small_primes() -> None:
    assert list(small_primes) == sorted(small_primes)

    assert not 0 in small_primes
    assert not 1 in small_primes
    assert 2 in small_primes

    for i in range(2, small_prime_limit):
        assert (i in small_primes) == all(i % f != 0 for f in range(2, i))

def test_modest_primes() -> None:
    for i in range(small_prime_limit):
        assert (i in modest_primes) == (i in small_primes)

    for i in range(small_prime_limit, min(small_prime_limit * small_prime_limit, modest_prime_limit)):
        assert (i in modest_primes) == all(i % p != 0 for p in small_primes)

def test_jacobi() -> None:
    import random
    for p in modest_primes:
        if p == 2:
            continue
        n = random.randint(1, p-1)
        pp = pow(n, p - 1 >> 1, p)
        if pp == p - 1:
            pp = -1
        assert jacobi(n, p) == pp

def test_is_square() -> None:
    for i in range(1000):
        assert is_square(i*i)
        for j in range(i*i + 1, (i+1)*(i+1)):
            assert not is_square(j)

def test_floor_sqrt() -> None:
    for i in range(1000):
        for j in range(i*i, i*i + 2*i + 1):
            assert floor_sqrt(j) == (i, i*i)

def test_floor_cbrt() -> None:
    for i in range(100):
        for j in range(i*i*i, i*i*i + 3*i*i + 3*i + 1):
            assert floor_cbrt(j) == (i, i*i*i)

def test_euclid() -> None:
    for i in range(1, 500):
        for j in range(1, 500):
            x, y, g = euclid(i, j)
            assert i % g == 0
            assert j % g == 0
            assert g == i * x + j * y
        x, y, g = euclid(i, 0)
        assert g == i
        assert x == 1
        x, y, g = euclid(0, i)
        assert g == i
        assert y == 1

def timecall(f: Any, *args: Any, **kwargs: Any) -> Tuple[float, Any]:
    import time
    start = time.time()
    r = f(*args, **kwargs)
    elapsed = time.time() - start
    return elapsed, r

if __name__ == '__main__':
    for offset in (-i * 0.0001 for i in range(14670,14680)):
        smp = sorted(p for p in modest_primes)
        upper = max(i/(p/(p.bit_length()+offset)) for i, p in enumerate(smp) if p > 500)
        lower = min(i/(p/(p.bit_length()+offset)) for i, p in enumerate(smp) if p > 500)
        print(offset, upper - lower, lower, upper)

    import time
    st = time.time()
    list(sieve_primes_to(1000000))
    print(time.time() - st)

    import cmath
    def rc(c: complex) -> complex:
        s = 1048576
        return round(c.real*s)/s + round(c.imag*s)/s*1j
    for p in range(3,500,2):
        if p == 2:
            continue
        S = 0j
        for k in range(1, p):
            S += jacobi(k, p) * cmath.rect(1, 2 * k * cmath.pi / p)
        print(p, rc(S*S), rc(S / sqrt(p)))
