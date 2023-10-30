from misc import tiny_primes, small_primes, small_primes_list

from math import gcd
from typing import Generator

def pollard_rho(N: int) -> int:
    # No attempt to behave sensibly if N is 1 or prime!
    for p in small_primes:
        if N % p == 0:
            return p
    # print(f'pollard_rho {N}')
    for i in range(5, N, 2):
        if N % i == 0:
            return i                    # Just in case!
        slow = i
        fast = i
        while True:
            slow = (slow * slow) % N + 1
            fast = (fast * fast) % N + 1
            fast = (fast * fast) % N + 1
            if slow == fast:
                break
            g = gcd(slow - fast, N)
            if g != 1:
                return g
    assert False, f'Failed on {N}'

def factor_repeats(N: int, factors: list[int]) -> Generator[int, None, None]:
    for f in factors:
        while N % f == 0:
            yield f
            N //= f
    assert N == 1

# TODO - this does the prime cert, that's excessive.
def is_prime(N: int) -> bool:
    if N in small_primes:
        return True
    for x in tiny_primes:
        if N % x == 0:
            return False
    if N < 2:
        return False
    reduced = N - 1
    twos = 0
    while (reduced & 1) == 0:
        reduced //= 2
        twos += 1
    # N is odd.
    factors = None
    for x in small_primes_list:
        if x * x > N:
            return True              # If we get here we have found any factors.
        # Note that if x is a factor of N then x^(N-1) is not 1 mod N.  Check
        # the squaring chain from x^reduced.  To certify we want x^((N-1)/2) ==
        # -1.
        p = pow(x, reduced, N)
        if p == 1:
            continue                    # Useless.
        useless = False
        for i in range(twos - 1):
            if p == N - 1:
                useless = True
                break                   # Useless, not a generator.
            q = p
            p = p * p % N
            if p == 1:
                # print(f'Extra sqrt {q} for {N}')
                return False            # Extra sqrt(1)
        if useless:
            continue
        if p != N - 1:
            #if p * p % N == 1:
                # print(f'Extra sqrt {p} for {N}')
            return False           # Either an extra sqrt(1) or incorrect order.

        # Now check the cofactors.
        if factors is None:
            factors = list(unique_prime_factors(reduced))
        for f in factors:
            if pow(x, (N - 1) // f, N) == 1:
                break                   # Not a generator.
        else:
            factors.sort()
            #print(f'{N} is prime gen {x} N-1 factors', ' '.join(str(f) for f in factors))
            return True
    assert False, 'Oops, failed to certify {N}'

def unique_prime_factors(N: int, verbose=True) -> Generator[int, None, None]:
    remain = N
    while remain > 1:
        if is_prime(remain):
            yield remain
            return
        f = pollard_rho(remain)
        g = remain // f
        if g < f:
            factor = g
            remain = f
        else:
            remain = g
            factor = f
        for p in unique_prime_factors(factor, False):
            #print(f'{N} has factor {p}')
            yield p
            while remain % p == 0:
                remain //= p
