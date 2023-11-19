
import constants, misc, pseudo_prime, pratt_cert

import array, joblib, math, mmap, sys
from typing import Iterator, Sequence, Tuple

sys.set_int_max_str_digits(1000000)

class FoundPrime(Exception): pass

def pocklington_serial(n: int) -> int:
    if n in constants.prime_by_bits:
        return constants.prime_by_bits[n]

    # The prime will be in the form p·t + 1 where t < p.  So first generate a
    # prime a bit bigger than 2^(n/2).
    p = pocklington_serial((n + 3) // 2)
    assert p * p >= 1 << n

    # Start with the lowest possible value of t.  We want t even, so ensure
    # that.  The value here is too small, the t+=2 below gets us in range.
    t = (1 << n-1) // p & -2

    while True:
        t += 2
        N = t * p + 1
        assert N.bit_length() == n, f'{n} {p} {t} {1<<n-1} {N} {1<<n}'
        if any(N != q and N % q == 0 for q in misc.small_primes):
            continue
        if pocklington_try_one(N, t, p):
            return N

def pocklington_parallel(n: int) -> int:
    if n in constants.prime_by_bits:
        return constants.prime_by_bits[n]
    if n < 1000:
        return pocklington_serial(n)

    # Generate an n bit prime (between 2^(n-1) and 2^n.  It will be in the form
    # p·t + 1 where t < p.  So generate a prime enough bigger than 2^(½n-½)
    # to give us some leeway.
    p = pocklington_parallel((n + 3) // 2)
    assert p*p >= 1 << n
    print(f'Search {n} bits', file=sys.stderr, flush=True)
    try:
        joblib.Parallel(n_jobs=-1, batch_size=1, timeout=86400)(
            joblib.delayed(pocklington_raise_one)(N, t, p)
            for N, t in pocklington_generate(p, n))
    except FoundPrime as e:
        return e.args[0]
    assert False                        # Hopefully we never get here!

def pocklington_generate(p: int, n: int) -> Iterator[Tuple[int, int]]:
    start_t = (1 << n-1) // p & -2
    # If we were smarter, we would sieve.
    plist = misc.small_primes if n < 2000 else misc.modest_primes_list
    # The big list takes ≈ 1 hour for all the reductions, so only do this when a
    # multi-hour runtime is expected.
    if n > 30000:
        try:
            plist = get_bertha()[0:100000000]     # Limit to 31 bits.
            print('Loading big bertha')
        except FileNotFoundError:
            pass
    # Pre-reduce p and smart_t modulo each element of plist.  Python doesn't
    # seem to give a nice way to parallelize this...
    p_mod = array.array('I')
    p_mod.extend(p % q for q in plist)
    start_t_mod = array.array('I')
    start_t_mod.extend(start_t % q for q in plist)
    print('Reduction done')
    # The range bound is arbitrary and should be quite a lot more than what is
    # necessary.
    for i in range(2, n * n * n, 2):
        t = start_t + i
        N = t * p + 1
        assert N.bit_length() == n, f'{n} {p} {t} {1<<n-1} {N} {1<<n}'

        # The condition is equivalent to `all(N % q != 0 for q in plist)` but
        # faster to compute.  Note that we never get N = q, as we only get
        # called with N ≥ 2^(n-1) ≫ q
        if all((pp * (ss + i) + 1) % q != 0
               for q, pp, ss in zip(plist, p_mod, start_t_mod)):
            yield N, t
        #if all(N % q != 0 for q in plist):
        #    yield N, t

def pocklington_raise_one(N: int, t: int, p: int) -> None:
    print(t & 65535, file=sys.stderr, flush=True)
    if pocklington_try_one(N, t, p):
        raise FoundPrime(N)

def pocklington_try_one(N: int, t: int, p: int) -> bool:
    # Pocklington test for N = p·t + 1 with N < p², p prime.  We believe the
    # caller on the primeness of p and computing N.
    assert p * p > N
    assert 1 <= t < p

    # Do a Fermat test for N.  It doesn't really matter what the base is, so
    # just use 2.
    if pow(2, N - 1, N) != 1:
        return False

    print('Passes fermat...', t & 65535, file=sys.stderr, flush=True)
    # If N is composite, then it has a prime factor q ≤ √N < p.
    #
    # As we pass the Fermat test 2^(N-1) ≡ 1 mod N, then 2^(N-1) ≡ 1 mod q.
    #
    # As q < p we have gcd(q-1, p) = 1, so that p is invertible mod (q-1), say
    # p·u ≡ 1 mod (q-1).
    #
    # Therefore (mod q), 2^t = 2^{(N-1)/p} ≡ 2^{(N-1)u} ≡ 1^u = 1.
    #
    # Now q divides both N and 2^t - 1, and so also the gcd below.
    #
    # Conversely, prime N will pass this test, unless 2^t ≡ 1 mod N.  (The
    # exceptional case should be rare: if N is prime, then there are t solutions
    # of x^t ≡ 1 mod N, and a randomly choosen x has probability 1/p of being an
    # exception).
    if math.gcd(pow(2, t, N) - 1, N) != 1:
        return False

    return True

bertha: Sequence[int]|None = None

def get_bertha() -> Sequence[int]:
    global bertha
    if bertha is not None:
        return bertha
    with open('primes.bin') as f:
        m = mmap.mmap(f.fileno(), 0, mmap.MAP_PRIVATE, mmap.ACCESS_READ)
        bertha = memoryview(m).cast('I')
    assert bertha[0] == 2
    assert bertha[-1] == 4294967291
    return bertha

# Assumes that the file exists.
def test_bertha() -> None:
    b = get_bertha()
    assert get_bertha() is b
    assert b[0] == 2
    assert b[-1] > 4000000000

if __name__ == '__main__':
    import time
    if len(sys.argv) > 1:
        for s in sys.argv[1:]:
            start = time.time()
            print(pocklington_parallel(int(s)))
            print('  took', time.time() - start, file=sys.stderr)
    else:
        def test_one(i: int) -> None:
            print(i)
            N = pocklington_parallel(i)
            print(N)
            assert pseudo_prime.baillie_psw(N)
        for i in range(2, 400):
            test_one(i)
        for i in 1000, 1001, 1002, 1500, 2000, 2500, 3000, 4000:
            test_one(i)
