
from misc import small_primes, small_primes_list, tiny_primes
from pollard_rho import pollard_rho, unique_prime_factors

from dataclasses import dataclass
from typing import Generator, Iterator, Optional, NamedTuple

@dataclass
class PrattCert:
    N: int
    is_prime: bool
    generator: int
    cofactors: list['PrattCert']

    def verify_no_rec(self):
        assert self.is_prime
        assert self.N > 1
        if False and self.N in small_primes:
            assert len(cofactors) == 0
            return

        assert 0 < self.generator < self.N
        order = self.N - 1
        remain = order
        for c in self.cofactors:
            assert remain % c.N == 0
            remain //= c.N
            while remain % c.N == 0:
                remain //= c.N
            assert pow(self.generator, order // c.N, self.N) != 1
        twos = 0
        while remain % 2 == 0:
            twos += 1
            remain //= 2
        assert remain == 1
        assert twos != 0
        assert pow(self.generator, (self.N - 1) // 2, self.N) == self.N - 1

    def verify(self):
        for c in self.cofactors:
            c.verify()
        self.verify_no_rec()

# Assumes that the iterator covers the sub-certs!
def verify_certs(d: dict[int, PrattCert]):
    for k, v in d.items():
        for c in v.cofactors:
            assert c.N in d
        v.verify_no_rec()

def pratt_cert(N: int, cache: dict[int,PrattCert]) -> Optional[PrattCert]:
    if N in cache:
        return cache[N]
    if N == 2:
        ret = PrattCert(N, True, 1, [])
        cache[N] = ret
        return ret
    if N == 1:
        return None
    twos = 0
    reduced = N - 1
    while (reduced & 1) == 0:
        reduced //= 2
        twos += 1
    if twos == 0:
        return None                     # N-1 is odd, N is even.
    cofactors = None
    for x in small_primes_list:
        x = x % N
        if x == 0:
            continue                    # Prime but we want a cert.
        p = pow(x, reduced, N)
        if p == 1:
            continue                    # Useless
        useless = False
        for i in range(twos - 1):
            if p == N - 1:
                useless = True
                break                   # Useless, not a generator.
            q = p
            p = p * p % N
            if p == 1:
                return None             # Extra sqrt(1).
        if useless:
            continue
        if p != N - 1:
            return None               # Either an extra sqrt(1) or fails Fermat.
        if cofactors is None:
            cofactors = list(unique_prime_factor_certs(reduced, cache))
        for f in cofactors:
            if pow(x, (N - 1) // f.N, N) == 1:
                break                   #  Not a generator.
        else:
            cofactors.sort(key=lambda c: c.N)
            print(f'{N} is prime gen {x} N-1 factors', ' '.join(
                str(f.N) for f in cofactors))
            ret = PrattCert(N, True, x, cofactors)
            cache[N] = ret
            return ret

    assert False, 'Oops, failed to certify {N}'

def unique_prime_factor_certs(N: int, cache: dict[int, PrattCert]) -> Generator[PrattCert, None, None]:
    remain = N
    while remain > 1:
        cert = pratt_cert(remain, cache)
        if cert:
            yield cert
            return
        f = pollard_rho(remain)
        g = remain // f
        if g < f:
            factor = g
            remain = f
        else:
            remain = g
            factor = f
        for p in unique_prime_factor_certs(factor, cache):
            yield p
            while remain % p.N == 0:
                remain //= p.N
