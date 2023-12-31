
import misc
from pollard_rho import pollard_rho, unique_prime_factors

from dataclasses import dataclass
from typing import Iterator, Optional, Tuple

@dataclass
class PrattCert:
    N: int
    generator: int
    cofactors: list['PrattCert']

    def verify_no_rec(self) -> None:
        if self.N == 2:
            assert self.generator == 1
            assert len(self.cofactors) == 0
            return
        assert self.N > 1
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

    def verify(self) -> None:
        for c in self.cofactors:
            c.verify()
        self.verify_no_rec()

# Assumes that the iterator covers the sub-certs!
def verify_certs(d: dict[int, PrattCert]) -> None:
    for k, v in d.items():
        for c in v.cofactors:
            assert c.N in d
        v.verify_no_rec()

def pratt_cert(N: int, cache: dict[int,PrattCert]) -> Optional[PrattCert]:
    if N in cache:
        return cache[N]
    if N == 2:
        ret = PrattCert(N, 1, [])
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
    for x in misc.small_primes:
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
        if N - p != 1:
            return None               # Either an extra sqrt(1) or fails Fermat.
        if cofactors is None:
            cofactors = list(unique_prime_factor_certs(reduced, cache))
        for f in cofactors:
            if pow(x, (N - 1) // f.N, N) == 1:
                break                   #  Not a generator.
        else:
            cofactors.sort(key=lambda c: c.N)
            # print(f'{N} is prime gen {x} N-1 factors', ' '.join(str(f.N) for f in cofactors))
            ret = PrattCert(N, x, cofactors)
            cache[N] = ret
            return ret

    assert False, 'Oops, failed to certify {N}'

# Uses pollard rho.
def unique_prime_factor_certs(N: int, cache: dict[int, PrattCert]) -> Iterator[PrattCert]:
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

def test_certs() -> None:
    certs: dict[int, PrattCert] = {}
    for p in misc.modest_primes:
        c = pratt_cert(p, certs)
        assert c
        assert c == certs[p]
    verify_certs(certs)

def test_composite() -> None:
    certs: dict[int, PrattCert] = {}
    for n in range(misc.modest_prime_limit):
        if not n in misc.modest_primes:
            assert pratt_cert(n, certs) is None
    verify_certs(certs)

if __name__ == '__main__':
    import sys
    certs: dict[int, PrattCert] = {}
    for s in sys.argv[1:]:
        print(f'{s}:', ' '.join(
            str(c.N) for c in unique_prime_factor_certs(eval(s), certs)))
    for prime in sorted(certs.keys()):
        cert = certs[prime]
        print(prime, 'gen', cert.generator, 'co',
              ' '.join(str(c.N) for c in cert.cofactors))
    verify_certs(certs)
