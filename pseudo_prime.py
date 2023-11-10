
import misc, quad

import math
from typing import Tuple

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

def miller_rabin1(n: int, base: int) -> bool:
    g = math.gcd(n, base)
    if g != 1:
        return g == n
    if n & 1 == 0:
        return n == 2
    Q, S = split_twos(n - 1)
    sq = pow(base, Q, n)
    if sq == 1:
        return True
    for i in range(S - 1):
        if n - sq == 1:
            return True                 # Got to -1.
        sq = sq * sq % n
        if sq == 1:                     # Unexpected √1,
            return False
    return n - sq == 1         # Otherwise either fails Fermat or unexpected √1.

# n must be odd, positive, and D=P²-4Q (mod n) with (D|n) = -1.
def strong_frobenius(n: int, P: int, Q: int, D: int) -> bool:
    #print(f'Frobenius {n=} {P=} {Q=} {D=} {n-D=}')
    P = P % n
    Q = Q % n
    # Don't reduce D: if it is small and negative, keep it that way.
    assert n > 0 and n & 1 == 1
    assert D % n == (P*P - 4*Q) % n
    assert misc.jacobi(D, n) == -1
    d, s = split_twos(n + 1)

    if P & 1:
        halfP = (P + n) >> 1
    else:
        halfP = P >> 1

    # ½P + ½√D
    base = quad.QuadInt(halfP, (1 + n) >> 1, quad.QuadRing(n, D))
    # Start from base^d, and follow through the chain of squares to base^(n+1).
    # Check that the square root of a rational is either rational or pure
    # quadratic.
    #
    # Then check that we get base^d = Q
    power = base.pow(d)

    if power.q == 0:
        # base^d is already a pure rational.  Only thing left to check is
        # base^(n+1).
        return pow(power.r, 1 << s, n) == Q

    # We don't have a pure rational.  Squaring enough times should get us there,
    # the step before should give a pure quadratic.
    twos = 0
    while power.r != 0:
        twos += 1
        if twos >= s:
            return False                # Went off into ga-ga land.
        power = power.square()
        if power.q == 0:
            # We found a rational square of a number that is neither pure
            # rational nor pure quadratic.
            return False

    # We got a pure quadratic.  Now finish the chain of squaring and check that
    # we end with Q.
    rational = power.q * power.q % n * D % n
    twos += 1

    # Now do the remainder of the squaring.
    return pow(rational, 1 << (s - twos), n) == Q

def strong_frobenius_a_star(n: int) -> bool:
    n = abs(n)
    if n & 1 == 0:
        return n == 2
    # Reject perfect squares: for a perfect square the check (D|n) = -1 never
    # passes and we would just do a linear scan looking for a factor!
    if n < 2 or misc.is_square(n):
        return False
    # Note that if the test below passes for a·b then it passes for one of a or
    # b.  However, we skip 3, so we cover multiples of 3 as well as primes.
    i = 5
    while True:
        if n == i:
            return True                 # We hit this for 5 and 11.
        if i & 3 == 3:
            D = -i
        else:
            D = i
        assert D % 4 == 1
        j = misc.jacobi(D, n)
        if j == 0:
            # i (or maybe 3) is a factor.  We never hit this on odd prime n.
            # Either we hit the n==i check above, or else n==3 and we choose
            # D=5.
            return False
        if j == -1:
            break
        i += 2

    if D == 5:                          # Special case to avoid Q = ±1.
        P = 5
        Q = 5
    else:
        P = 1
        Q = (1 - D) // 4

    return strong_frobenius(n, P, Q, D)

def baillie_psw(n: int) -> bool:
    if n < misc.modest_prime_limit:
        return n in misc.modest_primes
    if n.bit_length() < 1024:
        if any(n % p == 0 for p in misc.small_primes):
            return False
    else:
        if any(n % p == 0 for p in misc.modest_primes_list):
            return False
    return miller_rabin1(n, 2) and strong_frobenius_a_star(n)

def test_frob() -> None:
    # We need A* not just A for 5777...
    assert 5777 == 53 * 109
    assert strong_frobenius(5777, 1, -1, 5)
    assert not strong_frobenius_a_star(5777)

    assert 100981997 == 677 * 149161    # Ditto needs A*.
    assert strong_frobenius(100981997, 1, -1, 5)
    assert not strong_frobenius_a_star(100981997)

    assert strong_frobenius_a_star(1000981)
    assert strong_frobenius_a_star(65537)
    assert not strong_frobenius_a_star(5273095699)

def test_small() -> None:
    for n in range(misc.modest_prime_limit):
        assert strong_frobenius_a_star(n) == (n in misc.modest_primes)

def test_big():
    import constants
    assert baillie_psw(constants.prime_by_bits[512])
    assert baillie_psw(constants.prime_by_bits[2048])

def test_vpsp() -> None:
    # Pseudos (vpsp)!
    for n in 913, 150267335403, 430558874533, 14760229232131, 936916995253453:
        assert not strong_frobenius_a_star(n)
        if n == 14760229232131:
            P, Q, D = 1, 2, -7
        else:
            P, Q, D = 5, 5, 5
        assert misc.jacobi(D, n) == -1
        assert D == P*P - 4*Q
        ring = quad.QuadRing(n, D)
        base = quad.QuadInt(P, 1, ring) / quad.QuadInt(2, 0, ring)
        power = base.pow(n + 1)
        assert power.r == Q % n
        assert power.q != 0

if __name__ == '__main__':
    import sys, time
    sys.set_int_max_str_digits(1000000)
    st = time.time()
    for n in sys.argv[1:]:
        print(strong_frobenius_a_star(int(n)))
    print(time.time() - st)
