
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
        assert order & 1 == 0, f'{order} {self.N}'
        assert pratt_check_gen(self.N, self.cofactors, self.generator)
        assert pow(self.generator, order // 2, self.N) == self.N - 1

    def verify(self) -> None:
        for c in self.cofactors:
            c.verify()
        self.verify_no_rec()

# Use divide-and-conquer for computing the cofactor powers for the generator
# check.  This is the recursive function.
def pratt_check_worker(N: int, plist: list[Tuple[int, int]],
                       others: int, lo: int, hi: int) -> bool:
    if lo + 1 >= hi:
        p, e = plist[lo]
        return e <= 1 or pow(others, pow(p, e-1), N) != 1
    exp_lo = 1
    exp_hi = 1
    mid = (lo + hi) // 2
    for i in range(lo, mid):
        exp_lo *= pow(plist[i][0], plist[i][1])
    for i in range(mid, hi):
        exp_hi *= pow(plist[i][0], plist[i][1])
    # Note the swap!
    others_lo = pow(others, exp_hi, N)
    others_hi = pow(others, exp_lo, N)
    if others_lo == 1 or others_hi == 1: # Short circuit.
        return False
    return pratt_check_worker(N, plist, others_lo, lo, mid) \
        and pratt_check_worker(N, plist, others_hi, mid, hi)

# Check the generator for a Pratt cert.  The divide-and-conquer speeds the
# runtime for large numbers of cofactors.
def pratt_check_gen(N: int, cofactors: list[PrattCert], g: int) -> bool:
    if not cofactors:
        return True                     # Nothing to do.
    remain = N - 1
    plist = []
    for c in cofactors:
        e = 0
        while remain % c.N == 0:
            e += 1
            remain //= c.N
        assert e != 0
        plist.append((c.N, e))
    assert remain.bit_count() == 1
    assert remain > 1
    others = pow(g, remain, N)
    if others == 1:
        return False
    return pratt_check_worker(N, plist, others, 0, len(plist))

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
        if pratt_check_gen(N, cofactors, x):
            cofactors.sort(key=lambda c: c.N)
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
