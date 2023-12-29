
from dataclasses import dataclass, field
from typing import Tuple

@dataclass
class Monic:
    p: int
    # ∑ c_k x^k = 0, degree = len(c) -1 and c[degree] = 1.
    c: list[int]
    degree: int = 0
    small_pow: list[list[int]] = field(default_factory=list)

    def __post_init__(self) -> None:
        assert self.c
        assert self.c[-1] == 1
        self.degree = len(self.c) - 1
        for i in range(self.degree):
            self.small_pow.append([0] * i + [1])

        for _ in range(self.degree, 2 * self.degree):
            self.small_pow.append(self.multx(self.small_pow[-1]))

    def multx(self, ll: list[int]) -> list[int]:
        if len(ll) < self.degree:
            return [0] + ll

        h = ll[-1]
        r = [0] + ll[:-1]
        for i in range(self.degree):
            r[i] = (r[i] - h * self.c[i]) % self.p
        return r

    def multx_in_place(self, ll: list[int]) -> None:
        assert len(ll) == self.degree
        h = ll[-1]
        for i in reversed(range(1, self.degree)):
            ll[i] = (ll[i-1] - self.c[i] * h) % self.p
        ll[0] = self.c[0] * h % self.p

    def add(self, xx: list[int], yy: list[int]) -> list[int]:
        if len(xx) < len(yy):
            xx, yy = yy, xx
        r = list(xx)
        for i, y in enumerate(yy):
            r[i] = (r[i] + y) % self.p
        return r

    def sub(self, xx: list[int], yy: list[int]) -> list[int]:
        r = list(xx)
        if len(r) < len(yy):
            r.extend((0 for _ in range(len(yy) - len(r))))
        for i, y in enumerate(yy):
            r[i] = (r[i] - y) % self.p
        return r

    def mac_in_place(self, xx: list[int], v: int, yy: list[int]) -> None:
        assert len(yy) <= len(xx)
        for i, y in enumerate(yy):
            xx[i] = (xx[i] + v * y) % self.p

    def mult(self, xx: list[int], yy: list[int]) -> list[int]:
        # Multiply without reducing, then reduce.
        if xx == [] or yy == []:
            return []
        rr = [0] * (len(xx) + len(yy) - 1)
        #assert len(rr) <= len(self.small_pow), f'{self.degree} {len(xx)} {len(yy)} {len(self.small_pow)}'
        for i, x in enumerate(xx):
            for j, y in enumerate(yy):
                rr[i + j] += x * y
        if len(rr) <= self.degree:
            for i, r in enumerate(rr):
                rr[i] = r % self.p
            return rr

        result = rr[:self.degree]
        for i in range(self.degree, len(rr)):
            self.mac_in_place(result, rr[i], self.small_pow[i])

        return result

    def reduce(self, xx: list[int]) -> list[int]:
        r = xx[:self.degree]
        if len(r) < self.degree:
            r += [0] * (self.degree - len(r))
            return r
        for exp, v in enumerate(r[self.degree:self.degree*2], self.degree):
            self.mac_in_place(r, v, self.small_pow[exp])
        power = self.small_pow[-1]
        for v in r[self.degree*2:]:
            self.multx_in_place(power)
            self.mac_in_place(r, v, power)
        return r

    def pow(self, xx: list[int], n: int) -> list[int]:
        #print(f'Pow {self} {n}')
        if n == 0:
            return [1]
        result = xx
        for i in reversed(range(n.bit_length() - 1)):
            result = self.mult(result, result)
            if n & (1 << i):
                result = self.mult(result, xx)
        return result

    def powx(self, n: int) -> list[int]:
        if n == 0:
            return [1]
        result = [0, 1]
        for i in reversed(range(n.bit_length() - 1)):
            result = self.mult(result, result)
            if n & (1 << i):
                self.multx_in_place(result)
        return result

    def crack(self) -> list[int] | int | None:
        if self.degree <= 1:
            return None

        derivative = self.c[1:]
        for i, x in enumerate(derivative):
            derivative[i] = x * (i+1) % self.p
        _, _, gcd = euclid(self.c, derivative, self.p)
        if len(gcd) > 1:
            return gcd

        # Try for gcd with X^{p^i} - X for 1 ≤ 2i ≤ degree.
        X = self.small_pow[1]
        power = X                       # X¹ = X^p⁰
        # For even degree we need to try degree/2, for odd degree we need to
        # try (degree-1)/2.
        for i in range(1, (self.degree + 2) >> 1):
            power = self.pow(power, self.p) # X^{2^i}
            poly = list(power)
            # FIXME poly == 0 is possible here, if we split into polynomials
            # all of degree i.
            poly[1] = (poly[1] - 1) % self.p
            _, _, gcd = euclid(self.c, poly, self.p)
            if len(gcd) != 1:
                if len(gcd) < self.degree:
                    return gcd
                else:
                    assert self.degree % i == 0
                    return i
        return None

def poly_add(xx: list[int], yy: list[int], p: int) -> list[int]:
    if len(xx) < len(yy):
        xx, yy = yy, xx
    xx = list(xx)
    for i, y in enumerate(yy):
        xx[i] = (xx[i] + y) % p
    return xx

def poly_sub(xx: list[int], yy: list[int], p: int) -> list[int]:
    if len(xx) >= len(yy):
        xx = list(xx)
        for i, y in enumerate(yy):
            xx[i] = (xx[i] - y) % p
        return xx
    else:
        yy = list(yy)
        for i, x in enumerate(xx):
            yy[i] = (x - yy[i]) % p
        for i in range(len(xx), len(yy)):
            yy[i] = -yy[i] % p
        return yy

def poly_mult(xx: list[int], yy: list[int], p: int) -> list[int]:
    if not xx or not yy:
        return []
    r = [0] * (len(xx) + len(yy) - 1)
    for i, x in enumerate(xx):
        for j, y in enumerate(yy):
            k = i + j
            r[k] = (r[k] + x * y) % p
    return r

def poly_mult_scalar_in_place(xx: list[int], v: int, p: int) -> None:
    for i, x in enumerate(xx):
        xx[i] = x * v % p

def poly_div_mod(rr: list[int], yy: list[int], p: int) -> list[int]:
    # rr is reduced in place, quotient is returned.
    degree = len(yy) - 1
    if degree < 0:
        raise ZeroDivisionError('polynomial division by zero')
    while yy[degree] == 0:
        if degree == 0:
            raise ZeroDivisionError('polynomial division by zero')
        degree -= 1
    if len(rr) <= degree:
        return []
    monica = pow(yy[degree], -1, p)      # Possible div-by-zero!
    result = [0] * (len(rr) - degree)
    for i in reversed(range(len(result))):
        f = rr[i + degree] * monica % p
        result[i] = f
        rr[i + degree] = 0
        for j in range(degree):
            rr[i + j] = (rr[i + j] - f * yy[j]) % p
    #del rr[degree:]
    return result

def trim(xx: list[int]) -> None:
    if xx == []:
        return

    for degree in reversed(range(len(xx))):
        if xx[degree]:
            del xx[degree+1:]
            return

    del xx[:]

def euclid(xx: list[int], yy: list[int], p: int) -> Tuple[
        list[int], list[int], list[int]]:
    trim(xx)
    trim(yy)
    aa = list(xx)
    bb = list(yy)
    k_a_x = [1]
    k_a_y: list[int] = []
    k_b_x: list[int] = []
    k_b_y = [1]
    while aa:
        # prev_bb = list(bb)
        q = poly_div_mod(bb, aa, p)
        trim(bb)
        # assert prev_bb == trim(poly_add(poly_mult(q, aa, p), bb, p))
        k_b_x = poly_sub(k_b_x, poly_mult(q, k_a_x, p), p)
        k_b_y = poly_sub(k_b_y, poly_mult(q, k_a_y, p), p)
        bb, aa = aa, bb
        k_a_x, k_b_x = k_b_x, k_a_x
        k_a_y, k_b_y = k_b_y, k_a_y

    monique = pow(bb[-1], -1, p)
    poly_mult_scalar_in_place(k_b_x, monique, p)
    poly_mult_scalar_in_place(k_b_y, monique, p)
    poly_mult_scalar_in_place(bb, monique, p)
    return k_b_x, k_b_y, bb

def test_add_sub() -> None:
    a = [1,1,1,1,1,1]
    b = [1,1,1]

    assert poly_sub(a, b, 3) == [0,0,0,1,1,1]
    assert poly_sub(b, a, 3) == [0,0,0,2,2,2]

    print(poly_mult([1,1,1], [1,-1,1], 5))

def test_basic() -> None:
    m = Monic(11, [1, 0, 1, 0, 1])
    for i in range(1, len(m.small_pow)):
        assert m.small_pow[i] == m.multx(m.small_pow[i-1])
    print(f'{m.mult([1,1,1],[1,-1,1])=}')
    assert m.mult([1,1,1],[1,-1,1]) == [0,0,0,0]
    assert m.crack() == 2              # Splits into 'unknown' degree 2 factors.

def test_7() -> None:
    assert Monic(7, [4,3,4,1]).crack() is None
    assert Monic(7, [6,6,6,4,1]).crack() is None
    print(poly_mult([4,3,4,1], [6,6,6,4,1], 7))
    m = Monic(7, [3,0,3,1,4,4,1,1])
    assert m.mult([4,3,4,1], [6,6,6,4,1]) == [0] * 7
    assert m.crack() in ([4,3,4,1], [6,6,6,4,1])
