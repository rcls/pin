
import misc
import pseudo_prime
import pratt_cert

import math
import joblib
from dataclasses import dataclass
from typing import Iterator, Tuple
import sys

sys.set_int_max_str_digits(1000000)

class FoundPrime(Exception): pass

def pocklington_serial(n: int) -> int:
    # Generate an n bit prime.  I.e., 2^{n-1} ≤ p < 2^n.
    if n < 32:
        assert n > 1
        N = (1 << n-1) + 1
        cache: dict[int, pratt_cert.PrattCert] = {}
        while not pratt_cert.pratt_cert(N, cache):
            N += 2
        assert 1 << n-1 <= N < 1 << n
        return N

    # The prime will be in the form p·t + 1 where t < p.  So first generate a
    # prime a bit bigger than 2^(n/2).
    p = pocklington_serial((n + 3) // 2)
    assert p * p >= 1 << n

    # Start with the lowest possible value of t.  We want t even, so ensure
    # that.
    t = ((1 << n-1) // p) & -2

    while True:
        t += 2
        N = t * p + 1
        assert 1 << n-1 <= N < 1 << n, f'{n} {p} {t} {1<<n-1} {N} {1<<n}'
        assert p * p > N
        if pocklington_try_one(N, t, p):
            return N

def pocklington_parallel(n: int, runner: joblib.Parallel = None) -> int:
    if n < 1000:
        return pocklington_serial(n)

    if runner is None:
        runner = joblib.Parallel(n_jobs=-1, batch_size=1)

    # Generate an n bit prime (between 2^(n-1) and 2^n.  It will be in the from
    # p.t + 1 where t < p.  So generate a prime enough bigger than 2^((n-1)/2)
    # to give us some leeway.
    p = pocklington_parallel((n + 3) // 2, runner)
    assert p*p >= 1 << n
    print(f'Search {n} bits', file=sys.stderr, flush=True)
    try:
        runner(joblib.delayed(pocklington_raise_one)(N, t, p)
               for N, t in pocklington_generate(p, n))
    except FoundPrime as e:
        return e.args[0]
    assert False

def pocklington_generate(p: int, n: int) -> Iterator[Tuple[int, int]]:
    # If we were smarter, we would sieve.
    plist = misc.small_primes if n < 4096 else misc.modest_primes_list
    # We want t even.
    start_t = ((1 << n-1) // p) & -2
    for i in range(2, n * n, 2):
        t = start_t + i
        N = t * p + 1
        assert 1 << n-1 <= N < 1 << n, f'{n} {p} {t} {1<<n-1} {N} {1<<n}'
        assert p * p > N

        if all(N % q != 0 for q in plist):
            yield N, t

def pocklington_raise_one(N: int, t: int, p: int):
    print(t & 65535, file=sys.stderr, flush=True)
    if pocklington_try_one(N, t, p):
        raise FoundPrime(N)

def pocklington_try_one(N: int, t: int, p: int) -> bool:
    # Pocklington test for N = p·t + 1 with N < p², p prime.  Check t.  We
    # believe the caller on the primeness of p.
    assert 1 <= t < p

    # Do a Fermat test for N.  It doesn't really matter what the base is, so
    # just use 2.
    if pow(2, N - 1, N) != 1:
        return False

    print('Passes fermat...', t & 65535, file=sys.stderr, flush=True)
    # If N is composite, then it has a prime factor q ≤ √N < p.
    #
    # We have passed the Fermat test, 2^(q-1) ≡ 1 mod q.
    #
    # Also, gcd(q-1, p) = 1 as q < p, so that p is invertible mod (q-1), say p.u
    # ≡ 1 mod (q-1).
    #
    # Therefore (mod q), 2^t = 2^{(N-1)/p} ≡ 2^{(N-1)u} ≡ 1 as q|N and 2^{N-1} ≡
    # 1 (mod N).
    #
    # Now q divides both N and 2^t - 1, and so also the gcd below.
    #
    # Conversely, prime N will pass this test, unless 2^t ≡ 1 mod N.  (The
    # exceptional case should be rare: if N is prime, then there are at most t
    # solutions of x^t ≡ 1 mod N, so a randomly choosen x has probability ≤ 1/p
    # of being an exception).
    if math.gcd(pow(2, t, N) - 1, N) != 1:
        return False

    return True

if __name__ == '__main__':
    import sys, time
    if len(sys.argv) > 1:
        for s in sys.argv[1:]:
            start = time.time()
            print(pocklington_parallel(int(s)))
            print('  took', time.time() - start)
    else:
        def test_iterator(I):
            for i in I:
                print(i)
                N = pocklington_parallel(i)
                print(N)
                assert pseudo_prime.baillie_psw(N)
        test_iterator(range(2, 400))
        test_iterator([1000, 1001, 1002, 1500, 2000, 2500, 3000, 4000])

# Search 65536 bits started at t ≡ 31184 mod 65536.
