from misc import tiny_primes, small_primes
import pseudo_prime

from math import gcd
from typing import Iterator

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

def factor_repeats(N: int, factors: list[int]) -> Iterator[int]:
    for f in factors:
        while N % f == 0:
            yield f
            N //= f
    assert N == 1

def unique_prime_factors(N: int) -> Iterator[int]:
    remain = N
    while remain > 1:
        if pseudo_prime.baillie_psw(remain):
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
        for p in unique_prime_factors(factor):
            #print(f'{N} has factor {p}')
            yield p
            while remain % p == 0:
                remain //= p
