#!/usr/bin/python3

import misc
from ratlinstring import ratlinstring

from fractions import Fraction
from math import sqrt, gcd
from numbers import Rational
from dataclasses import dataclass

@dataclass(frozen=True)
class Quadratic:
    # r + q × √b
    r: Rational
    q: Rational = Fraction(0)
    b: int = 5
    def reduce(self) -> 'Quadratic':
        b = self.b
        q = self.q
        # In leiu of a proper factorisation...
        for p in misc.tiny_primes:
            pp = p * p
            while b % pp == 0:
                b //= pp
                q *= p
        return Quadratic(self.r, q, b)
    def __float__(self) -> float:
        return float(self.r) + float(self.q) * sqrt(self.b)
    def __add__(self, y: 'Quadratic') -> 'Quadratic':
        assert self.b == y.b
        return Quadratic(self.r + y.r, self.q + y.q, self.b)
    def __sub__(self, y: 'Quadratic') -> 'Quadratic':
        assert self.b == y.b
        return Quadratic(self.r - y.r, self.q - y.q, self.b)
    def __neg__(self) -> 'Quadratic':
        return Quadratic(-self.r, -self.q, self.b)
    def __mul__(self, y: 'Quadratic') -> 'Quadratic':
        assert self.b == y.b
        return Quadratic(self.r * y.r + self.q * y.q * self.b,
                         self.r * y.q + self.q * y.r, self.b)
    def __bool__(self) -> bool:
        return self.r != 0 or self.q != 0
    def invert(self) -> 'Quadratic':
        norm = (self.r * self.r) - self.q * self.q * self.b
        return Quadratic(self.r / norm, -self.q / norm, self.b)
    def __truediv__(self, y: 'Quadratic') -> 'Quadratic':
        return self * y.invert()
    def addi(self, n: int) -> 'Quadratic':
        return Quadratic(self.r + Fraction(n), self.q, self.b)
    def multi(self, n: int) -> 'Quadratic':
        return Quadratic(self.r * n, self.q * n, self.b)
    def __repr__(self) -> str:
        return ratlinstring(self.r, self.q, f'√{self.b} ')

# ax+b / cx+d
@dataclass(frozen=True)
class Moïbus:
    a: int
    b: int
    c: int
    d: int
    def apply(self, x: float) -> float:
        return (self.a * x + self.b) / (self.c * x + self.d)
    def applyi(self, n: int) -> Fraction:
        return Fraction(self.a * n + self.b, self.c * n + self.d)
    def applyf(self, x: Rational) -> Fraction:
        return Fraction(self.a * x.numerator + self.b * x.denominator,
                        self.c * x.numerator + self.d * x.denominator)
    def applyq(self, x: Quadratic) -> Quadratic:
        return (x.multi(self.a).addi(self.b)) / \
            (x.multi(self.c).addi(self.d))
    # We deal with determinant ±1 so no need to normalise...
    def __mul__(self, y: 'Moïbus') -> 'Moïbus':
        return Moïbus(self.a * y.a + self.b * y.c,
                      self.a * y.b + self.b * y.d,
                      self.c * y.a + self.d * y.c,
                      self.c * y.b + self.d * y.d)
    def fixed_point(self) -> Quadratic:
        # Solve x = self(x).
        # x = (a·x + b) / (cx + d)
        # c·x² + (d - a)·x - b = 0.
        # Take +ve sqrt.
        ad = self.a - self.d
        common = gcd(ad, self.b * 2, self.c * 2)
        base = (ad * ad + 4 * self.b * self.c) // (common * common)
        sqbase = int(sqrt(base))
        if base == sqbase * sqbase:
            # No need for a quad. term!
            result = Quadratic(
                Fraction(ad + common * sqbase, self.c * 2), Fraction(0), 0)
        else:
            result = Quadratic(Fraction(ad, self.c * 2),
                               Fraction(common, self.c * 2), base)
        assert result == self.applyq(result), \
            f'Fixed point {self} gives {result} to {self.applyq(result)}'
        return result

    def __repr__(self) -> str:
        return f'({self.a}𝑥 + {self.b})/({self.c}𝑥 + {self.d})'

def addconst(x: int) -> Moïbus:
    return Moïbus(1, x, 0, 1)

ident = Moïbus(1, 0, 0, 1)
recip = Moïbus(0, 1, 1, 0)
zero = Moïbus(0, 0, 0, 1)
inf = Moïbus(0, 1, 0, 0)

def chain(l: list[Moïbus]) -> Moïbus:
    r = ident
    for x in l:
        r = r * x
    return r

def approx(x: float, tolerance: float = 1e-7) -> Quadratic:
    f = ident
    coeffs: list[int] = []
    items: list[Moïbus] = []
    residues: list[float] = []
    r = x
    while True:
        residues.append(r)
        n = int(r)
        frac = r - n
        approx = chain(items).applyf(Fraction(n))
        if frac == 0 or abs(float(approx) - x) <= tolerance :
            return chain(items).applyq(Quadratic(Fraction(n), Fraction(0), 0))
        # Look for a quad. within tolerance.
        for i in range(len(items)):
            res = chain(items[:i]).applyq(chain(items[i:]).fixed_point())
            if abs(float(res) - x) <= tolerance:
                return res.reduce()
        r = 1 / frac
        coeffs.append(n)
        item = addconst(n) * recip
        items.append(item)
        f = f * item
        #print(f, x - f.apply(1), x - (f * inf).apply(0), r, coeffs)

if __name__ == '__main__':
    import sys
    import math
    from math import sqrt, exp, log, e, pi
    if len(sys.argv) > 1:
        if len(sys.argv) == 2:
            n = 2
            tol = 1e-7
        else:
            n = len(sys.argv) - 1
            tol = float(sys.argv[-1])
        for v in sys.argv[1:n]:
            x = float(eval(v)) # x = float(v)
            a = approx(x, tol)
            av = float(a)
            err = av - x
            print(f'{v} ≈ {a} [err = {err}]')
            if a.b == 5:
                r = a.r - a.q
                q = a.q * Fraction(2, 1)
                print('    =', ratlinstring(r, q, 'φ'))
    else:
        print(approx(sqrt(2)))
        print(approx(0.61803399, 5e-9))
        print(approx(pi, 1e-10))
        print(approx(e, 1e-10))
