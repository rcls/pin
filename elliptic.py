
import misc

from dataclasses import dataclass
import math
from typing import Iterator, Tuple

Point = Tuple[int, int]

@dataclass
class Elliptic:
    # y² ≡ x³ + a·x + b (mod p).
    a: int
    # b is never actually needed except when we verify a point is valid.
    b: int
    p: int
    negp: int                           # Python sucks.

    @staticmethod
    def curve(x, y, a, p) -> 'Elliptic':
        if a >= p or a <= p // -2:
            a = a % p
        assert a != 0
        b = (y * y - x * x % p * x - a * x) % p
        assert b != 0
        return Elliptic(a, b, p, -p)

    def verify(self, z: Point) -> None:
        x, y = z
        if x == 0 == y:
            return
        assert (x * x % self.p * x + self.a * x + self.b - y * y) % self.p == 0

    def neg(self, z: Point) -> Point:
        x, y = z
        return x, -y

    def canon(self, z: Point) -> Point:
        x, y = z
        return x % self.p, y % self.p

    def add(self, z0: Point, z1: Point) -> Point:
        x0, y0 = z0
        x1, y1 = z1
        if x0 == 0 == y0:
            return z1
        if x1 == 0 == y1:
            return z0
        # Parameterise as x = (1 - u)·x₀ + u·x₁, y = (1 - u)·y₀ + u·y₁
        # i.e., x = u·(x₁ - x₀) + x₀, y = u·(y₁ - y₀) + y₀,
        #
        # Express x³ - y² + a·x + b ≡ 0 in terms of u.
        xs = x1 - x0
        ys = y1 - y0
        if xs != 0 and xs != self.p and xs != self.negp:
            # Normal case.
            identical = False
            pass
        elif ys == 0 or ys == self.p or ys == self.negp:
            identical = True
            # In the exceptional case where z0 == z1, instead we take a tangent.
            #
            # Differentiating: 2·y·dy = (3·x² + a)·dx.  So take the slopes:
            xs = y0 * 2
            ys = (x0 * x0 * 3 + self.a) % self.p
        elif (y0 + y1) % self.p == 0:
            # The identity…  Represent as (0,0) which is not an otherwise valid
            # point.
            return 0, 0
        else:
            # x₀ ≡ x₁ and y₀ ≢ ±y₁.  We have y₀² ≡ y₁², and so a factorisation
            # of p.
            g = math.gcd(ys, self.p)
            assert 1 < g < self.p
            # Maybe one day we will hit this?
            print(self.p, g)
            assert False
            raise misc.FoundFactor(self.p, g)

        # x³  = xs³·u³ + 3·xs²·x₀·u² + 3·xs·x₀²·u + x₀³
        # y²  =               ys².u² + 2·ys·y₀ ·u + y₀²
        # a·x =                        a·xs    ·u + a·x0
        # b = b
        xs2 = xs * xs % self.p
        # Note that c2 is not reduced.
        c3 = xs2 * xs % self.p
        c2 = xs2 * (x0 * 3) - ys * ys
        #c1 = x0x0_3_p_a * xs - ys * (y0 * 2)
        if c3 == 0:
            # xs³ is a multiple of p, so we must have a common factor.
            g = math.gcd(xs, self.p)
            assert 1 < g < self.p
            raise misc.FoundFactor(self.p, g)
        # Special handling for the indentical case case again...
        if identical:
            # Now our polynomial has a double root at u=0 and the solution
            # we want us c₃·u + c₂ = 0
            d = c2
        else:
            # c₀ = x₀³ - y₀² + a·x₀ + b ≡ 0 by assumption.
            #
            # The polynomial equation is c₃·u³ + c₂·u² + c₁·u ≡ 0, which is
            # satisfied by u=0 and u=1.  (The latter implies, c₃ + c₂ + c₁ ≡ 0.)
            #assert (c1 + c2 + c3) % self.p == 0
            #
            # Dividing by u gives c₃·u² + c₂·u + c₁, and then dividing by u-1
            # gives c₃·u + c₃ + c₂, so u = -(c₃ + c₂) / c₃.
            d = c3 + c2

        try:
            uneg = d % self.p * pow(c3, -1, self.p) % self.p
        except ValueError:
            g = math.gcd(self.p, c3)
            assert 1 < g < self.p
            raise misc.FoundFactor(self.p, g)

        # Now compute x and y at u = -uneg, and flip the sign of y.
        return (x0 - xs * uneg) % self.p,  (ys * uneg - y0) % self.p

    def double(self, z: Point) -> Point:
        # Specialize add() for the doubling case.
        x, y = z
        if x == 0 == y:
            return z

        xs = y * 2
        ys = (x * x * 3 + self.a) % self.p

        xs2 = xs * xs % self.p
        c3 = xs2 * xs % self.p
        c2 = xs2 * (x * 3) - ys * ys

        if c3 == 0:
            # xs³ is a multiple of p, so we must have a common factor.
            g = math.gcd(xs, self.p)
            assert 1 < g < self.p, f'{xs}, {self}'
            raise misc.FoundFactor(self.p, g, xs)

        try:
            uneg = c2 % self.p * pow(c3, -1, self.p) % self.p
        except ValueError:
            g = math.gcd(self.p, c3)
            assert 1 < g < self.p
            raise misc.FoundFactor(self.p, g)
        return (x - xs * uneg) % self.p,  (ys * uneg - y) % self.p

    def mult(self, z: Point, n: int) -> Point:
        if n == 0 or z == (0, 0):
            return 0, 0
        if n < 0:
            n = -n
            z = (z[0], -z[1])
        result = (0,0)
        power2 = z
        # Do the calculation bottom up.  For the ECM, n (mod 2ᵏ) passes through
        # all possible values reasonably often, giving us more chance of hitting
        # repeated factors of the order.
        for i in range(n.bit_length() - 1):
            if n & (1 << i):
                result = self.add(result, power2)
            power2 = self.double(power2)
        return self.add(power2, result)

def test_basic() -> None:
    curve = Elliptic.curve(4, 5, 7, 65537)
    z = (4, 5)
    curve.verify(z)
    print(curve)
    z2 = curve.add(z, z)
    curve.verify(z2)
    z3 = curve.add(z, z2)
    curve.verify(z3)
    z4 = curve.add(z, z3)
    curve.verify(z4)
    assert curve.add(z2, z) == z3
    assert curve.add(z3, z) == z4
    assert curve.add(z2, z2) == z4
    nz = curve.neg(z)
    curve.verify(nz)
    assert curve.add(z2, nz) == z
    assert curve.add(z3, nz) == z2
    assert curve.add(z4, nz) == z3
    assert curve.add(z3, curve.neg(z4)) == curve.canon(nz)

def test_identity() -> None:
    import pytest
    curve = Elliptic.curve(4, 5, 7, 65537)
    z = (4, 5)
    nz = curve.neg(z)
    curve.verify(nz)
    assert curve.add(z, nz) == (0, 0)
    assert curve.add(z, (0,0)) == z
    assert curve.add((0,0), z) == z

def test_mult() -> None:
    curve = Elliptic.curve(4, 5, 7, 65537)
    z = (4, 5)
    direct = (0, 0)
    for i in range(1, 22):
        direct = curve.add(direct, z)
        curve.verify(direct)
        assert direct == curve.mult(z, i), f'{i}'

def test_find_order() -> None:
    curve = Elliptic.curve(4, 5, 7, 65537)
    z = (4, 5)
    # The order should be between 65537 ± 512
    zb = curve.mult(z, 65537 - 512)
    for i in range(1024):
        if zb == (0, 0):
            print('Order is', 65537 - 512 + i)
        zb = curve.add(z, zb)

def qtry_one(n: int, B1: int) -> None:
    import random
    a = random.randint(1,1000000) % n
    x = random.randint(1,1000000) % n
    y = random.randint(1,1000000) % n
    curve = Elliptic.curve(x, y, a, n)
    Q = (x,y)
    print('[', flush=True, end='')
    for p in misc.sieve_primes_to(B1):
        pp = p
        while pp < B1:
            Q = curve.mult(Q, p)
            pp *= p
    B2 = B1 * 100
    B1 -= B1 % 30                       # Adjust to a multiple of 30.
    if B1 == 0:
        return
    print(']', flush=True, end='')
    try:
        prod = 1
        Q1 = Q
        Q2 = curve.double(Q1)
        Q4 = curve.double(Q2)
        Q7 = curve.add(Q4, curve.add(Q2, Q1))
        Q11 = curve.add(Q7, Q4)
        Q13 = curve.add(Q11, Q2)
        Q30 = curve.add(Q13, curve.add(Q13, Q4))

        Q = curve.mult(Q30, B1 // 30)
        for _ in range(B1, B2, 6):
            # Todo - we could center on 15 mod 30 and use 2,4,8,14.
            prod = prod * (Q[0] - Q1[0]) % curve.p
            prod = prod * (Q[0] - Q7[0]) % curve.p
            prod = prod * (Q[0] - Q11[0]) % curve.p
            prod = prod * (Q[0] - Q13[0]) % curve.p
            Q = curve.add(Q, Q30)
    except:
        print('!')
        raise
    # TODO - detect if we get a zero above somewhere.
    g = math.gcd(prod, curve.p)
    assert g != curve.p
    if 1 < g < curve.p:
        print('*')
        raise misc.FoundFactor(curve.p, g)

def qtry_ret(n: int, l: int) -> int:
    try:
        qtry_one(n, l)
        return 0
    except misc.FoundFactor as e:
        print(f'Found factor {e.args[1]} at {l=}')
        return e.args[1]

_SCHEDULE = (
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

def qtry_single(n: int) -> int:
    for b1 in schedule():
        f = qtry_ret(n, b1)
    print('Found factor', f, 'of', n)
    return f

def qtry_parallel(n: int) -> None:
    import joblib
    try:
        joblib.Parallel(n_jobs=-1, batch_size=1)(
            joblib.delayed(qtry_one)(n, b1) for b1 in schedule())
    except misc.FoundFactor as e:
        print('Found factor', e.args[1], 'of', e.args[0], flush=True)

def qtry_parallel10(n: int) -> None:
    import joblib
    jobs = joblib.Parallel(n_jobs=-1, batch_size=1, return_as='generator')(
        joblib.delayed(qtry_ret)(n, l) for l in schedule())
    count = 0
    for f in jobs:
        if f != 0:
            count +=1
            if count == 10:
                break

def retry(curve, z, i):
    try:
        m = curve.mult(z, i)
        print('  retry succeeds')
        return m
    except misc.FoundFactor:
        print('  retry fails')
        return None

if __name__ == '__main__':
    #qtry_parallel((1<<101)-1)
    #qtry_parallel((1<<137)-1)
    qtry_parallel((1<<149)-1)
