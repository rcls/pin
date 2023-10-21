from misc import small_primes
from typing import Optional

class Quad:
    p: int                              # Should be odd prime
    S: int
    Q: int
    _non_residue: Optional[int]
    _square_chain: Optional[list[int]] = []

    def __init__(self, p, non_residue=None):
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
        if self.p == 2 or n % self.p == 0:
            return True
        assert self.p & 1 == 1
        r = pow(n, self.p >> 1, self.p)
        assert r == 1 or r == self.p - 1
        return r == 1

    def non_residue(self) -> int:
        p = self.p
        assert p & 1 == 1
        if self._non_residue:
            return self._non_residue
        for nr in small_primes:
            if nr == p:
                continue                # Useless.
            tp = pow(nr, p >> 1, p)
            if p - tp == 1:
                self._non_residue = nr
                return nr
            assert tp == 1, f'{p} is not prime (Fermat fail)'

        assert False, f'Failed to find non-residue for {p}'

    def square_chain(self, squares: int) -> int:
        if not self._square_chain:
            self._square_chain = [pow(self.non_residue(), self.Q, self.p)]
        while len(self._square_chain) <= squares:
            x = self._square_chain[-1]
            self._square_chain.append(x * x % self.p)
        return self._square_chain[squares]

    def sqrt(self, n: int) -> Optional[int]:
        p = self.p
        n = n % p
        if n == 0:
            return 0                    # Special case!
        if not self.is_qr(n):
            return None

        assert self.Q & 1 == 1

        root = pow(n, (self.Q + 1) >> 1, p)  # Square root of something.
        overshoot = pow(n, -1, p) * root % p * root % p
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

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        bs = sys.argv[-1]
        b = Quad(eval(bs))
        for s in sys.argv[1:-1]:
            v = eval(s)
            print(f'sqrt({s}) mod {bs} = {b.sqrt(v)}')
    else:
        q137 = Quad(137)
        print(f'{q137.is_qr(1)=}')
        print(f'{q137.sqrt(1)=}')
        print(f'{q137.sqrt(136)=}')
        print(f'{q137.is_qr(2)=}')
        print(f'{q137.sqrt(2)=}')
        print(f'{q137.is_qr(3)=}')
        print(f'{Quad(139).is_qr(-1)=}')
        q12s64p1 = Quad((12 << 64) + 1)
        print(f'{q12s64p1.sqrt(2)=}')
        print(f'{q12s64p1.sqrt(3)=}')
        print(f'{q12s64p1.sqrt(-3)=}')
        print(f'{q12s64p1.sqrt(5)=}')
        print(f'{q12s64p1.sqrt(-1)=}')

        q3s64m1 = Quad((3 << 64) - 1)
        print(f'{q3s64m1.sqrt(3)=}')
        print(f'{q3s64m1.sqrt(-3)=}')
