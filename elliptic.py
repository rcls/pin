#!/usr/bin/python3

# On 4096 bit numbers, the curve operations are about 6 times faster than
# Weierstrauss.  But at 150 bits, I think the python overheads dominate.

import misc

import math
from dataclasses import dataclass
from typing import Iterator, Tuple

# Projective coordinates (X:Z) with x = X/Z.
@dataclass(frozen=True)
class Point:
    X: int
    Z: int

# Raised if our computations go wrong without a factorization.
class Unexpected(Exception): pass

class MontgomeryCurve:
    p: int
    A: int
    Am2d4: int                          # (A-2)/4
    # We never use B, or the y co-ordinates.  For validation we only care
    # about whether or not B·y² ≡ x³ + A·x² + x is a q.r. or not.  j is
    # the Jacobi symbol (B|p).  [For composite p, we use the Jacobi symbol
    # rather than strictly whether or not B is a q.r.]
    j: int

    def invmod(self, x):
        return pow(x, -1, self.p)

    def moddiv(self, x, y):
        # x/y mod p
        return x * pow(y, -1, self.p) % self.p

    def __init__(self, p: int, Am2d4: int, j: int):
        self.p = p
        self.Am2d4 = Am2d4
        self.A = Am2d4 * 4 + 2
        self.gcd_check(Am2d4)
        self.gcd_check(Am2d4 + 1)
        assert self.A % p != 2 and -self.A % p != 2
        self.j = j

    @staticmethod
    def curve_on(p: int, Am2d4: int, x: int, z: int = 1):
        A = Am2d4 * 4 + 2
        j = misc.jacobi((x * x % p + A * x + 1) * x, p) * misc.jacobi(z, p)
        return MontgomeryCurve(p, Am2d4, j)

    def a_validate(self, x: int) -> None:
        p = self.p
        assert misc.jacobi((x * x % p + self.A * x + 1) * x, p) == self.j

    def validate(self, P: Point) -> None:
        # If we cared we could avoid the modular inversion by doing a projective
        # calculation.
        self.a_validate(self.reduce(P))

    # Extract a factorization after something has gone wrong arithemetically.
    def gcd_check(self, *args: int) -> None:
        for x in args:
            g = math.gcd(x, self.p)
            if 1 < g < self.p:
                raise misc.FoundFactor(self.p, g)

    # Map the projective coordinates to affine.
    def reduce(self, P: Point) -> int:
        return self.moddiv(P.X, P.Z)

    def eq(self, P: Point, Q: Point) -> int:
        return (P.X * Q.Z - P.Z * Q.X) % self.p == 0

    def a_add(self, a: int, b: int, n: int) -> int:
        # s·n = (a·b - 1)² / (a - b)².
        # s = (a·b - 1)² / [n·(a - b)²]
        p = self.p
        abm1 = a * b % p - 1
        amb = a - b
        return self.moddiv(abm1 * abm1 % p, amb * amb % p * n % p)

    def a_double(self, x:int) -> int:
        # We intentionally do little to simplify:
        #
        # (x² - 1)² / (4(x³ + Ax² + x))
        p = self.p

        xm1 = x * x - 1
        num = xm1 * xm1 % p

        den = ((x + self.A) * x % p + 1) * x * 4 % p

        return self.moddiv(num, den)

    # Add P+Q given N = P-Q.
    def add(self, P: Point, Q: Point, N: Point) -> Point:
        # Special cases.  The first two are essential to apply if we hit them,
        # the others are optional.
        if N.Z == 0:
            return self.double(P)       # P-Q = O, P=Q, P+Q = 2P
        if N.X == 0:
            D = self.double(P)
            return Point(D.Z, D.X)      # P-Q = T, P+Q = T + 2P = T + 2Q
        if P.Z == 0:
            return Q                    # P = O, P+Q = Q
        if Q.Z == 0:
            return P                    # Q = O, P+Q = P

        if P.X == 0:
            return Point(Q.Z, Q.X)      # P=T, invert Q.
        if Q.X == 0:
            return Point(P.Z, P.X)      # Q=T, invert P.

        X0, Z0 = P.X, P.Z
        X1, Z1 = Q.X, Q.Z

        # The affine coordinate of the result r is given by
        #
        # r·n = (p·q - 1)² / (p - q)².

        # Doing this projectively, we get four multiplies and two squarings:
        #
        # r·n = (X₀·X₁ - Z₀·Z₁)² / (X₀·Z₁ - Z₀·X₁)²
        #
        # But we can save some multiplications:
        #
        # 2(p·q - 1) = (p - 1)(q + 1) + (p + 1)(q - 1)
        # 2(p   - q) = (p - 1)(q + 1) - (p + 1)(q - 1)
        #
        # Noting the shared addends, the RHS above has just two multiplies, even
        # in projective form, replacing four in our first version.  The factors
        # of 2 cancel and so can be ignored.
        p = self.p
        PmQp = (X0 - Z0) * (X1 + Z1) % p
        PpQm = (X0 + Z0) * (X1 - Z1) % p

        num = PmQp + PpQm
        den = PmQp - PpQm

        r = Point(num * num % p * N.Z % p, den * den % p * N.X % p)
        if r.X == 0 and num != 0 and num != p:
            self.gcd_check(num, N.Z)
        if r.Z == 0 and den != 0:
            self.gcd_check(den, N.X)

        if r.X == 0 and r.Z == 0:
            raise Unexpected(N, P, Q)
        return r

    def double(self, P: Point) -> Point:
        # Special cases are P = O or T.  In both cases the result is O.
        if P.X == 0 or P.Z == 0:
            return Point(1, 0)

        # In affine coordinates: (x² - 1)² / (4(x³ + Ax² + x))
        #
        # The numerator is (x+1)²·(x-1)².  We can reuse (x+1)² to rewrite the
        # denominator as:
        #
        #     4x · [(x + 1)² + ¼(A-2)·4x]
        #
        # To homogenize, we need to replace 4x with 4XZ.  Again, reusing the
        # factors of the numerator, we can use the equality 4x = (x+1)² - (x-1)²
        # to avoid the explicit multiplication X·Z.
        p = self.p
        X, Z = P.X, P.Z
        XpZ = X + Z
        XmZ = X - Z
        XpZ2 = XpZ * XpZ % p
        XmZ2 = XmZ * XmZ % p
        fourXZ = XpZ2 - XmZ2

        # We assume ¼(A-2) is small, and so skip the reduction of ¼(A-2)·4XZ.
        r = Point(XpZ2 * XmZ2 % p,
                  fourXZ * (XpZ2 + self.Am2d4 * fourXZ) % p)
        if r.X != 0 or r.Z != 0:
            return r

        # Check everything...
        print('Collide')
        self.gcd_check(X, Z, XpZ, XmZ)
        raise Unexpected(P)               # Something wierd happened.

    def ladder2(self, n: int, P: Point) -> Tuple[Point, Point]:
        # Return ([n]P, [n+1]P).
        if n == 0:
            return Point(1, 0), P

        R, S = P, self.double(P)
        #r, s = 1, 2

        # Note that we skip the most significant 1 bit.
        for i in reversed(range(n.bit_length() - 1)):
            if n & (1 << i) == 0:
                R, S = self.double(R), self.add(R, S, P)
                #r, s = 2*r, r+s
            else:
                R, S = self.add(R, S, P), self.double(S)
                #r, s = r+s, 2*s

        #assert r == n and s == r + 1
        return R, S

    def ladder_mult(self, n: int, P: Point) -> Point:
        if n == 0:
            return Point(1, 0)

        # Short circuit factors of 2, we're not doing crypto.
        while n & 1 == 0:
            n >>= 1
            P = self.double(P)

        if n == 1:
            return P

        R, S = self.ladder2(n // 2, P)
        # If n=2m+1 is odd, n = m+(m+1)
        return self.add(R, S, P)

def test_basic() -> None:
    curve = MontgomeryCurve(65537, 25, -1)

    P = Point(3, 5)
    x = curve.reduce(P)
    curve.a_validate(x)

    P2 = curve.double(P)
    x2 = curve.a_double(x)
    assert curve.reduce(P2) == x2
    curve.a_validate(x2)

    P3 = curve.add(P, P2, P)
    x3 = curve.a_add(x, x2, x)
    assert curve.reduce(P3) == x3
    curve.a_validate(x3)

    P4 = curve.double(P2)
    x4 = curve.a_double(x2)
    assert curve.reduce(P4) == x4
    curve.a_validate(x4)

    P4a = curve.add(P, P3, P2)
    x4a = curve.a_add(x, x3, x2)
    assert curve.reduce(P4a) == x4a
    assert x4a == x4

def test_ladder() -> None:
    curve = MontgomeryCurve(65537, 27, -1)
    O = Point(1, 0)
    P = Point(3, 4)
    b = curve.reduce(P)
    multx = [O, P, curve.double(P)]
    multa = [0, b, curve.a_double(b)]
    for i in range(3, 1000):
        multx.append(curve.  add(multx[-1], P, multx[-2]))
        multa.append(curve.a_add(multa[-1], b, multa[-2]))

    for i in range(1, 1000):
        assert curve.reduce(multx[i]) == multa[i]

    for i in range(2, 1000, 2):
        assert curve.eq(multx[i], curve.double(multx[i//2]))
        assert multa[i] == curve.a_double(multa[i//2])

    for i in range(1, 1000):
        L = curve.ladder_mult(i, P)
        assert multa[i] == curve.reduce(L), f'{i}'

def test_order() -> None:
    curve = MontgomeryCurve(65537, 27, -1)
    P = Point(3, 4)
    base = 65536 - 512
    Q0, Q1 = curve.ladder2(base, P)
    assert curve.eq(Q0, curve.ladder_mult(base, P))
    order = None
    for i in range(0, 2048):
        Q0, Q1 = Q1, curve.add(Q1, P, Q0)
        if Q1.Z == 0:
            order = base + i + 2
            print(i, order, Q1, curve.ladder_mult(order, P))
    assert order is not None

def test_order_long() -> None:
    curve = MontgomeryCurve(65537, 27, -1)
    P = Point(3, 4)
    Q0, Q1 = Point(1, 0), P
    order = None
    for i in range(2, 65536 + 512):
        Q0, Q1 = Q1, curve.add(Q1, P, Q0)
        if Q1.Z == 0:
            order = i
            print(i, Q1, curve.ladder_mult(i, P))
    assert order is not None

def btry_one(n: int, B1: int, Am2d4: int):
    import time
    print('.', end='', flush=True)
    start1 = time.time()
    #import random
    #Am2d4 = random.randint(1, n - 2)
    Q = Point(2,1)                      # Doesn't really matter?
    curve = MontgomeryCurve.curve_on(n, Am2d4, Q.X, Q.Z)
    for p in misc.sieve_primes_to(B1):
        pp = p
        ppp = pp * p
        while ppp < B1:
            pp = ppp
            ppp *= p
        Q = curve.ladder_mult(pp, Q)

    curve.gcd_check(Q.X, Q.Z)

    if Q.X == 0 or Q.Z == 0:
        print('Trivial!')
        return

    start2 = time.time()
    try:
        Q2 = curve.double(Q)
        Q3 = curve.add(Q2, Q, Q)
        Q4 = curve.add(Q3, Q, Q2)
        Q6 = curve.add(Q4, Q2, Q2)
        Q7 = curve.add(Q4, Q3, Q)
        Q11 = curve.add(Q7, Q4, Q3)
        Q13 = curve.add(Q7, Q6, Q)
        Q15 = curve.add(Q13, Q2, Q11)
        Q30 = curve.double(Q15)

        B2 = B1 * 100
        R0, R1 = curve.ladder2(B1//30, Q30)

        check = Q.Z * Q7.Z % curve.p * Q11.Z % curve.p * Q13.Z % curve.p
        if check == 0:
            curve.gcd_check(Q.Z, Q7.Z, Q11.Z, Q13.Z)
            raise Unexpected(Q.Z, Q7.Z, Q11.Z, Q13.Z)

        for _ in range(B1, B2, 30):
            R0, R1 = R1, curve.add(R1, Q30, R0)
            for L in Q, Q7, Q11, Q13:
                det = (L.X * R1.Z - L.Z * R1.X) % curve.p
                if det == 0:
                    curve.gcd_check(check)
                    print('Stage 2 trivial')
                    return
                # No need to include L.Z again!
                #
                # Do we actually need the R1.Z check here?  The chance of it
                # hitting appears negligable as R1 = [30n]Q & we should have
                # already covered the factors of 30n.
                ncheck = check * det % curve.p #* R1.Z % curve.p
                if ncheck == 0:
                    curve.gcd_check(check, R1.Z, det)
                    raise Unexpected(check, R1.Z, det)
                check = ncheck
    except:
        print('!', flush=True)
        raise
    try:
        curve.gcd_check(check)
    except:
        print('*', flush=True)
        raise
    end = time.time()
    #print(start2 - start1, end - start2)

_SCHEDULE = (
    #(5, 100, 4),
    #(10, 500, 10),
    (15,       2000,     25),
    (20,      11000,     90),
    (25,      50000,    300),
    (30,     250000,    700),
    (35,    1000000,   1800),
    (40,    3000000,   5100),
    (45,   11000000,  10600),
    (50,   43000000,  19300),
    (55,  110000000,  49000),
    (60,  260000000, 124000),
    (65,  850000000, 210000),
    (70, 2900000000, 340000))

def schedule(min_digits:int = 0) -> Iterator[int]:
    for digits, B1, count in _SCHEDULE:
        if digits >= min_digits:
            for _ in range(count):
                yield B1
    while True:                         # Not that we actually get here...
        yield B1

def btry_parallel(n: int, min_digits:int = 0) -> int:
    import joblib
    try:
        joblib.Parallel(n_jobs=-1, batch_size=1)(
            joblib.delayed(btry_one)(n, b1, sqAm2d4 * (sqAm2d4 + 1))
            for sqAm2d4, b1 in enumerate(schedule(min_digits), 1))
        assert False
    except misc.FoundFactor as e:
        print('Found factor', e.args[1]) #, 'of', e.args[0], flush=True)
        return e.args[1]

def drive(n: int) -> int|None:
    for p in misc.modest_primes_list:
        if p * p > n:
            return None
        if n % p == 0:
            return p
    if n < misc.modest_prime_limit * misc.modest_prime_limit:
        return None

    return btry_parallel(n)

if __name__ == '__main__':
    import sys, pseudo_prime
    sys.set_int_max_str_digits(1<<30)
    if len(sys.argv) == 2:
        S = sys.argv[1]
        N = eval(S)
        if pseudo_prime.baillie_psw(N):
            print('Probable prime')
        else:
            print('Factor', S)
            btry_parallel(N)
        sys.exit(0)
    elif len(sys.argv) == 3:
        S = sys.argv[1]
        N = eval(S)
        btry_parallel(N, eval(sys.argv[2]))
        sys.exit(0)
    import math, mersenne
    for p in misc.modest_primes_list:
        if p <= 1277:
            continue
        if p > 1500:
            break
        if p in (523, 727, 751, 809, 971, 983, 997, 1061, 1277, 1237):
            continue
        if mersenne.Mersenne(p).is_prime():
            continue
        N = (1<<p) - 1
        print(f'M({p}) = {N}')
        f = drive(N)
        assert f is not None
        print(f'M({p}) has factor {f} (bits: {math.log2(f):.2f})')
    # 1061:473
    # 1237:230
    # 1277
    # 1283:132

    # 727:323
    # 523:226
    # 751:217
    # 809:201
    # 997:187
    # 971:174
    # 983:140
    #
    # 293:85
    # 563:81
    # 647:78
    # 739:77
    # 599:74
    # 347:74
    # 149:66
    # 137:65
    # 433:64
    # (math.factorial(62) - 1)//13143173//458476324671361
    # ...352600227238301502049958597925403
    #btry_parallel(4378942815107578007)
    #btry_one((1<<149) - 1, 2000)
    #pass
