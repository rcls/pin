import misc, ratlinstring

from dataclasses import dataclass

class QuadRing:
    p: int                              # Should be odd prime
    S: int
    Q: int
    non_residue: int
    _square_chain: list[int] = []

    def __init__(self, p: int, non_residue: int|None = None):
        # TODO - how to fail if p<2.
        self.p = p
        self.Q, self.S = misc.split_twos(p - 1)
        if non_residue is not None:
            self.non_residue = non_residue
            return

        assert p & 1 == 1
        assert not misc.is_square(p)
        for nr in misc.modest_primes_list:
            if nr == p:
                continue
            j = misc.jacobi(nr, p)
            assert j != 0, f'{p} is not prime, factor {j}'
            if j == -1:
                self.non_residue = nr
                break
        else:
            assert False

    def is_qr(self, n: int) -> bool:
        # Only valid if p really is prime!
        return self.p == 2 or misc.jacobi(n, self.p) >= 0

    def square_chain(self, squares: int) -> int:
        if not self._square_chain:
            self._square_chain = [pow(self.non_residue, self.Q, self.p)]
        while len(self._square_chain) <= squares:
            x = self._square_chain[-1]
            self._square_chain.append(x * x % self.p)
        return self._square_chain[squares]

    def maybe_sqrt(self, n: int) -> int|None:
        p = self.p
        n = n % p
        if n == 0:
            return 0                    # Special case!
        if not self.is_qr(n):
            return None

        assert self.Q & 1 == 1

        root = pow(n, (self.Q + 1) >> 1, p)  # Square root of something.
        overshoot = pow(n, -1, p) * root % p * root % p # ≡ n^{Q+1} / n = n^Q.
        assert overshoot == pow(n, self.Q, p)
        twos = self.S - 1
        # From the fact that n is a q.r., we have that the multiplicative order
        # of overshoot is a power of two, bounded by self.S-1.
        assert pow(overshoot, 1 << twos, p) == 1
        assert root * root % p == overshoot * n % p
        while overshoot != 1:
            raised = overshoot
            # Find `twos` such that overshoot^(2^twos) ≡ -1.  This should get a
            # smaller on each successive iteration.
            for twos in range(twos):
                if p - raised == 1:
                    break
                raised = raised * raised % p
                assert raised != 1
            else:
                assert False, 'Failed to steer...'

            # Find multiplier such that multiplier^(2^(twos+1)) ≡ -1,
            # i.e., multiplier^(2^twos) is a square-root of -1.
            multiplier = self.square_chain(self.S - 2 - twos)
            assert p - pow(multiplier, 2 << twos, p) == 1
            # Multiply root by multiplier, and hence to maintain our invariants,
            # multiply overshoot by multiplier².  This reduces the
            # multiplicative order of overshoot.
            overshoot = overshoot * multiplier % p * multiplier % p
            root = root * multiplier % p
            assert pow(overshoot, 1 << twos, p) == 1
            assert root * root % p == overshoot * n % p
        return root

    # Use at own risk, no checking...
    def sqrt(self, n: int) -> int:
        p = self.p
        n = n % p
        if n == 0:
            return 0                    # Special case!

        root = pow(n, (self.Q + 1) >> 1, p)  # Square root of something.
        overshoot = pow(n, -1, p) * root % p * root % p # ≡ n^(Q+1)/n = n^Q.
        twos = self.S - 1
        while overshoot != 1:
            raised = overshoot
            # We should need a smaller value for twos on each iteration.
            for twos in range(twos):
                if p - raised == 1:
                    break
                raised = raised * raised % p
            else:
                assert False, 'Failed to steer...'

            multiplier = self.square_chain(self.S - 2 - twos)
            overshoot = overshoot * multiplier % p * multiplier % p
            root = root * multiplier % p
        return root

    def qsqrt(self, n: int) -> 'QuadInt':
        if self.is_qr(n):
            return QuadInt(self.sqrt(n), 0, self)
        else:
            p, k = self.p, self.non_residue
            return QuadInt(0, self.sqrt(n * k) * pow(k, -1, p) % p, self)

    def __str__(self) -> str:
        return f'√{self.non_residue} (mod {self.p})'

@dataclass(frozen=True)
class QuadInt:
    # r + q√k
    r: int
    q: int
    k: QuadRing
    def __init__(self, r:int, q:int, k: QuadRing):
        object.__setattr__(self, 'r', r % k.p)
        object.__setattr__(self, 'q', q % k.p)
        object.__setattr__(self, 'k', k)
    def __add__(self, y: 'QuadInt') -> 'QuadInt':
        assert self.k == y.k
        return QuadInt(self.r + y.r, self.q + y.q, self.k)
    def __sub__(self, y: 'QuadInt') -> 'QuadInt':
        assert self.k == y.k
        return QuadInt(self.r - y.r, self.q - y.q, self.k)
    def __neg__(self) -> 'QuadInt':
        return QuadInt(-self.r, -self.q, self.k)
    def __mul__(self, y: 'QuadInt') -> 'QuadInt':
        r, q, k = self.r, self.q, self.k
        assert k == y.k
        return QuadInt(r * y.r + q * y.q % k.p * k.non_residue,
                       r * y.q + q * y.r, k)
    def __truediv__(self, y: 'QuadInt') -> 'QuadInt':
        return self.__mul__(y.invert())
    def square(self) -> 'QuadInt':
        r, q, k = self.r, self.q, self.k
        D = k.non_residue
        if D.bit_length() < 32:
            return QuadInt(r * r + q * q * D, 2 * r * q, k)
        else:
            return QuadInt(r * r + q * q % k.p * D, 2 * r * q, k)

    def invert(self) -> 'QuadInt':
        # (r - q√k)(r² - q²k)⁻¹
        r, q, p = self.r, self.q, self.k.p
        modsq = r * r - q * q % p * self.k.non_residue
        invmodsq = pow(modsq, -1, p)
        return QuadInt(r * invmodsq, -q * invmodsq, self.k)

    def halve(self) -> 'QuadInt':
        r, q, p = self.r, self.q, self.k.p
        r = r // 2 if r & 1 == 0 else (r + p) // 2
        q = q // 2 if q & 1 == 0 else (q + p) // 2
        return QuadInt(r, q, self.k)

    def pow(self, n: int) -> 'QuadInt':
        #print(f'Pow {self} {n}')
        if n == 0:
            return QuadInt(1, 0, self.k)
        b = self
        if n < 0:
            b = b.invert()
            n = -n
        result = b
        for i in reversed(range(n.bit_length() - 1)):
            result = result.square()
            if n & (1 << i):
                result *= b
        return result

    def __str__(self) -> str:
        from fractions import Fraction
        return ratlinstring.ratlinstring(
            Fraction(self.r, 1), Fraction(self.q, 1),
            f' √{self.k.non_residue} ') + f' (mod {self.k.p})'

# Test helper.
def check_square_sqrt(q: QuadRing, x: int) -> None:
    s = q.maybe_sqrt(x)
    assert s is not None
    assert misc.jacobi(x, q.p) >= 0
    assert s == q.sqrt(x)
    ss = s * s % q.p
    assert ss == x % q.p
    #if s * 2 > q.p:
    #    s -= q.p

def test_137() -> None:
    q137 = QuadRing(137)
    print(f'{q137=!s}')
    assert q137.is_qr(1)
    assert q137.maybe_sqrt(1) == 1

    sqrtm1 = q137.maybe_sqrt(-1)
    assert sqrtm1 is not None
    check_square_sqrt(q137, sqrtm1)
    assert q137.is_qr(2)
    check_square_sqrt(q137, 2)
    assert misc.jacobi(3, 137) == -1
    assert not q137.is_qr(3)
    assert q137.maybe_sqrt(3) == None
    print(f'{q137=!s}')

    assert not QuadRing(139).is_qr(-1)

def test_12_64() -> None:
    q12s64p1 = QuadRing((12 << 64) + 1)
    print(f'{q12s64p1=}')
    check_square_sqrt(q12s64p1, 2)
    check_square_sqrt(q12s64p1, 3)
    check_square_sqrt(q12s64p1, -3)
    assert q12s64p1.maybe_sqrt(5) is None
    assert misc.jacobi(5, q12s64p1.p) == -1
    check_square_sqrt(q12s64p1, -1)

def test_3_64() -> None:
    q3s64m1 = QuadRing((3 << 64) - 1)
    assert misc.jacobi(2, q3s64m1.p)
    check_square_sqrt(q3s64m1, 2)
    print(f'{q3s64m1.qsqrt(2)=!s}')
    print(f'{q3s64m1.maybe_sqrt(3)=}')
    print(f'{q3s64m1.qsqrt(-3)=!s}')
    check_square_sqrt(q3s64m1, 3)
    def square(x: QuadInt) -> QuadInt: return x*x
    print(f'{-square(q3s64m1.qsqrt(-3))=!s}')
    assert -square(q3s64m1.qsqrt(-3)) == QuadInt(3, 0, q3s64m1)
    print(q3s64m1)
    r = q3s64m1.p * q3s64m1.p - 1
    print(f'{r:#x}', f'{QuadInt(3,4,q3s64m1).pow(r)=!r}')
    assert QuadInt(3,4,q3s64m1).pow(r) == QuadInt(1, 0, q3s64m1)

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        bs = sys.argv[-1]
        b = QuadRing(eval(bs))
        for s in sys.argv[1:-1]:
            v = eval(s)
            print(f'sqrt({s}) mod {bs} = {b.maybe_sqrt(v)}')
