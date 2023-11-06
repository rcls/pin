from misc import small_primes, jacobi

from dataclasses import dataclass
from typing import Optional

class Quad:
    p: int                              # Should be odd prime
    S: int
    Q: int
    _non_residue: Optional[int]
    _square_chain: Optional[list[int]] = []

    def __init__(self, p: int, non_residue: Optional[int] = None):
        # TODO - how to fail if p<2.
        self.p = p
        self._non_residue = non_residue
        Q = p - 1
        S = 0
        while Q & 1 == 0:
            Q >>= 1
            S += 1
        self.S = S
        self.Q = Q

    def is_qr(self, n: int) -> bool:
        # Only valid if p really is prime!
        return self.p == 2 or jacobi(n, self.p) >= 0

    def non_residue(self) -> int:
        p = self.p
        assert p & 1 == 1
        if self._non_residue:
            return self._non_residue
        for nr in small_primes:
            if nr == p:
                continue                # Useless.
            j = jacobi(nr, p)
            assert j != 0, f'{p} is not prime, factor {j}'
            if j == -1:
                self._non_residue = nr
                return nr

        assert False, f'Failed to find non-residue for {p}'

    def square_chain(self, squares: int) -> int:
        if not self._square_chain:
            self._square_chain = [pow(self.non_residue(), self.Q, self.p)]
        while len(self._square_chain) <= squares:
            x = self._square_chain[-1]
            self._square_chain.append(x * x % self.p)
        return self._square_chain[squares]

    def maybe_sqrt(self, n: int) -> Optional[int]:
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
        # From the fact that n is a q.r...
        assert pow(overshoot, 1 << twos, p) == 1
        assert root * root % p == overshoot * n % p
        while overshoot != 1:
            raised = overshoot
            # We should need a smaller value for twos on each iteration.
            for twos in range(twos):
                if p - raised == 1:
                    break
                raised = raised * raised % p
                assert raised != 1
            else:
                assert False, 'Failed to steer...'

            multiplier = self.square_chain(self.S - 2 - twos)
            assert p - pow(multiplier, 2 << twos, p) == 1
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
            p, k = self.p, self.non_residue()
            return QuadInt(0, self.sqrt(n * k) * pow(k, -1, p) % p, self)

    def __str__(self) -> str:
        if self._non_residue:
            return f'√{self._non_residue} (mod {self.p})'
        else:
            return f'√� (mod {self.p})'

@dataclass(frozen=True)
class QuadInt:
    # r + q√k
    r: int
    q: int
    k: Quad
    def __init__(self, r:int, q:int, k: Quad):
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
        return QuadInt(r * y.r + q * y.q % k.p * k.non_residue(),
                       r * y.q + q * y.r, k)
    def square(self) -> 'QuadInt':
        r, q, k = self.r, self.q, self.k
        D = k.non_residue()
        if D < 32:
            return QuadInt(r * r + q * q * D, 2 * r * q, k)
        else:
            return QuadInt(r * r + q * q % k.p * D, 2 * r * q, k)

    def invert(self) -> 'QuadInt':
        # (r - q√k)(r² - q²k)⁻¹
        r, q, p = self.r, self.q, self.k.p
        modsq = r * r - q * q % p * self.k.non_residue()
        invmodsq = pow(modsq, -1, p)
        return QuadInt(r * invmodsq, -q * invmodsq, self.k)

    def pow(self, n: int) -> 'QuadInt':
        #print(f'Pow {self} {n}')
        b = self
        if n == 0:
            return QuadInt(1, 0, self.k)
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
        return f'{self.r} + {self.q} {self.k}'

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        bs = sys.argv[-1]
        b = Quad(eval(bs))
        for s in sys.argv[1:-1]:
            v = eval(s)
            print(f'sqrt({s}) mod {bs} = {b.maybe_sqrt(v)}')
    else:
        q137 = Quad(137)
        print(f'{q137=!s}')
        assert q137.is_qr(1)
        assert q137.maybe_sqrt(1) == 1
        def square_sqrt(q, x):
            s = q.maybe_sqrt(x)
            assert s is not None
            assert jacobi(x, q.p) >= 0
            assert s == q.sqrt(x)
            ss = s * s % q.p
            assert ss == x % q.p
            if s * 2 > q.p:
                s -= q.p
            return s

        square_sqrt(q137, q137.maybe_sqrt(-1))
        assert q137.is_qr(2)
        assert square_sqrt(q137, 2)
        assert jacobi(3, 137) == -1
        assert not q137.is_qr(3)
        assert q137.maybe_sqrt(3) == None
        print(f'{q137=!s}')

        assert not Quad(139).is_qr(-1)

        q12s64p1 = Quad((12 << 64) + 1)
        square_sqrt(q12s64p1, 2)
        square_sqrt(q12s64p1, 3)
        square_sqrt(q12s64p1, -3)
        assert q12s64p1.maybe_sqrt(5) is None
        assert jacobi(5, q12s64p1.p) == -1
        square_sqrt(q12s64p1, -1)
        print(q12s64p1)

        q3s64m1 = Quad((3 << 64) - 1)
        assert jacobi(2, q3s64m1.p)
        square_sqrt(q3s64m1, 2)
        print(f'{q3s64m1.qsqrt(2)=!s}')
        print(f'{q3s64m1.maybe_sqrt(3)=}')
        def square(x): return x*x
        print(f'{q3s64m1.qsqrt(-3)=!s}')
        print(f'{-square(q3s64m1.qsqrt(-3))=!s}')
        print(q3s64m1)
        r = q3s64m1.p * q3s64m1.p - 1
        print(f'{r:#x}', f'{QuadInt(3,4,q3s64m1).pow(r)=!r}')
