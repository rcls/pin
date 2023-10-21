#!/usr/bin/python3

from math import sqrt, gcd
from dataclasses import dataclass

@dataclass(frozen=True)
class Fraction:
    n: int
    d: int
    def __init__(self, n: int, d: int = 1):
        g = gcd(n, d)
        if (d < 0):
            g = -g
        object.__setattr__(self, 'n', n // g)
        object.__setattr__(self, 'd', d // g)
    def value(self) -> float:
        return self.n / self.d
    def __mul__(self, y: 'Fraction') -> 'Fraction':
        return Fraction(self.n * y.n, self.d * y.d)
    def __truediv__(self, y: 'Fraction') -> 'Fraction':
        return Fraction(self.n * y.d, self.d * y.n)
    def __add__(self, y: 'Fraction') -> 'Fraction':
        return Fraction(self.n * y.d + self.d * y.n, self.d * y.d)
    def __sub__(self, y: 'Fraction') -> 'Fraction':
        return Fraction(self.n * y.d - self.d * y.n, self.d * y.d)
    def __neg__(self) -> 'Fraction':
        return Fraction(-self.n, self.d)
    def multi(self, n: int) -> 'Fraction':
        return Fraction(self.n * n, self.d)
    def __str__(self) -> str:
        if self.d == 1:
            return str(self.n)
        s = '-' if self.n < 0 else '';
        an = abs(self.n)
        if an == 1 and self.d <= 10:
            return s + '01½⅓¼⅕⅙⅐⅛⅑⅒'[self.d]
        if an + 1 == self.d and self.d <= 7:
            return '0½⅔¾⅘⅚'[an]
        if self.d == 5 and an < 5:
            return s + '0⅕⅖⅗⅘'[an]
        if self.d == 8 and an < 8:
            return s + '0⅛¼⅜½⅝¾⅞'[an]
        return f'{self.n}/{self.d}'

@dataclass(frozen=True)
class Quadratic:
    # r + q * sqrt(b)
    r: Fraction
    q: Fraction
    b: int
    def reduce(self) -> 'Quadratic':
        b = self.b
        q = self.q
        # In leiu of a proper factorisation...
        for p in 2, 3, 5, 7, 11, 13, 17, 19, 23, 29:
            pp = p * p
            while b % pp == 0:
                b //= pp
                q = q.multi(p)
        return Quadratic(self.r, q, b)
    def value(self) -> float:
        return self.r.value() + self.q.value() * sqrt(self.b)
    def multf(self, y: Fraction) -> 'Quadratic':
        ya = Fraction(abs(y.d), abs(y.n))
        return Quadratic(self.r * y, self.q * ya * y, self.b)
    def __add__(self, y: 'Quadratic') -> 'Quadratic':
        assert self.b == y.b
        return Quadratic(self.r + y.r, self.q + y.q, self.b)
    def __mul__(self, y: 'Quadratic') -> 'Quadratic':
        assert self.b == y.b
        return Quadratic(self.r * y.r + (self.q * y.q).multi(self.b),
                         self.r * y.q + self.q * y.r, self.b)
    def invert(self) -> 'Quadratic':
        norm = (self.r * self.r) - (self.q * self.q).multi(self.b)
        return Quadratic(self.r / norm, -self.q / norm, self.b)
    def __truediv__(self, y: 'Quadratic') -> 'Quadratic':
        return self * y.invert()
    def addi(self, n: int) -> 'Quadratic':
        return Quadratic(self.r + Fraction(n), self.q, self.b)
    def multi(self, n: int) -> 'Quadratic':
        return Quadratic(self.r.multi(n), self.q.multi(n), self.b)
    def __str__(self) -> str:
        if self.q.n == 0:
            return str(self.r)
        aq = Fraction(abs(self.q.n), abs(self.q.d))
        if aq.d == aq.n:
            straq =f'√{self.b}'
        elif aq.n == 1:
            straq = f'√{self.b} / {aq.d}'
        else:
            straq = f'{aq} × √{self.b}'
        if self.r.n == 0:
            return '-' + straq if self.q.n < 0 else straq
        s = '-' if self.q.n < 0 else '+'
        return f'{self.r} {s} {straq}'

# ax+b / cx+d
@dataclass(frozen=True)
class Rational:
    a: int
    b: int
    c: int
    d: int
    def apply(self, x: float) -> float:
        return (self.a * x + self.b) / (self.c * x + self.d)
    def applyi(self, n: int) -> Fraction:
        return Fraction(self.a * n + self.b, self.c * n + self.d)
    def applyf(self, x: Fraction) -> Fraction:
        return Fraction(self.a * x.n + self.b * x.d,
                        self.c * x.n + self.d * x.d)
    def applyq(self, x: Quadratic) -> Quadratic:
        return (x.multi(self.a).addi(self.b)) / \
            (x.multi(self.c).addi(self.d))
    # We deal with determinate ±1 so no need to normalise...
    def __mul__(self, y: 'Rational') -> 'Rational':
        return Rational(self.a * y.a + self.b * y.c,
                        self.a * y.b + self.b * y.d,
                        self.c * y.a + self.d * y.c,
                        self.c * y.b + self.d * y.d)
    def fixed_point(self) -> Quadratic:
        # x = self(x).
        # x = (ax + b) / (cx + d)
        # c x² + (d - a)x - b = 0.
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

def addconst(x: int) -> Rational:
    return Rational(1, x, 0, 1)

ident = Rational(1, 0, 0, 1)
recip = Rational(0, 1, 1, 0)
zero = Rational(0, 0, 0, 1)
inf = Rational(0, 1, 0, 0)

def chain(l: list[Rational]) -> Rational:
    r = ident
    for x in l:
        r = r * x
    return r

def approx(x: float, tolerance: float = 1e-7) -> Quadratic:
    f = ident
    coeffs: list[int] = []
    items: list[Rational] = []
    residues: list[float] = []
    r = x
    while True:
        residues.append(r)
        n = int(r)
        frac = r - n
        approx = chain(items).applyf(Fraction(n))
        if frac == 0 or abs(approx.value() - x) <= tolerance :
            return chain(items).applyq(Quadratic(Fraction(n), Fraction(0), 0))
        # Look for a quad. within tolerance.
        for i in range(len(items)):
            res = chain(items[:i]).applyq(chain(items[i:]).fixed_point())
            if abs(res.value() - x) <= tolerance:
                print('Quadratic', i, coeffs, res)
                return res.reduce()
        r = 1 / frac
        coeffs.append(n)
        item = addconst(n) * recip;
        items.append(item)
        f = f * item
        #print(f, x - f.apply(1), x - (f * inf).apply(0), r, coeffs)

if __name__ == '__main__':
    import sys
    from math import e, pi, sqrt
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
            av = a.value()
            err = av - x
            print(f'{v} ≈ {a} [err = {err}]')
            if a.b == 5:
                r = a.r - a.q
                q = a.q * Fraction(2, 1)
                print(f'    = {r} + {q} φ')
    else:
        print(approx(sqrt(2)))
        print(approx(0.61803399, 5e-9))
        print(approx(pi, 1e-10))
        print(approx(e, 1e-10))

