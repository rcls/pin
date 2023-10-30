
from dataclasses import dataclass, field
from typing import Union, Tuple

@dataclass
class Monic:
    p: int
    # x^order + ∑ c_k x^k = 0 where order = len(c)
    c: list[int]
    order: int = 0
    small_pow: list[list[int]] = field(default_factory=list)

    def __post_init__(self):
        self.order = len(self.c)
        for i in range(self.order):
            r = [0] * self.order
            r[i] = 1
            self.small_pow.append(r)
        for _ in range(self.order, 2 * self.order):
            self.small_pow.append(self.multx(self.small_pow[-1]))

    def multx(self, ll: list[int]) -> list[int]:
        h = ll[-1]
        r = [0] + ll[:-1]
        for i in range(self.order):
            r[i] = (r[i] - h * self.c[i]) % self.p
        return r
    def multx_in_place(self, ll: list[int]):
        h = ll[-1]
        for i in reversed(range(1, self.order)):
            ll[i] = (ll[i-1] - self.c[i] * h) % self.p
        ll[0] = self.c[i] * h % self.p

    def add(self, xx: list[int], yy: list[int]) -> list[int]:
        assert len(xx) == self.order
        assert len(yy) == self.order
        return [(x + y) % self.p for x, y in zip(xx, yy)]
    def sub(self, xx: list[int], yy: list[int]) -> list[int]:
        assert len(xx) == self.order
        assert len(yy) == self.order
        return [(x - y) % self.p for x, y in zip(xx, yy)]

    def add_in_place(self, xx: list[int], yy: list[int]):
        assert len(xx) == self.order
        assert len(yy) == self.order
        for i in range(len(yy)):
            xx[i] = (xx[i] + yy[i]) % self.p
    def mac_in_place(self, xx: list[int], v: int, yy: list[int]):
        assert len(xx) == self.order
        assert len(yy) == self.order
        for i in range(self.order):
            xx[i] = (xx[i] + v * yy[i]) % self.p

    def mult(self, xx: list[int], yy: list[int]) -> list[int]:
        # Multiply without reducing, then reduce.
        if xx == [] or yy == []:
            return []
        rr = [0] * (len(xx) + len(yy) - 1)
        #assert len(rr) <= len(self.small_pow), f'{self.order} {len(xx)} {len(yy)} {len(self.small_pow)}'
        for i, x in enumerate(xx):
            for j, y in enumerate(yy):
                rr[i + j] += x * y
        if len(rr) <= self.order:
            for i, r in enumerate(rr):
                rr[i] = r % self.p
            return rr

        reduce = rr[:self.order]
        for i in range(self.order, len(rr)):
            self.mac_in_place(reduce, rr[i], self.small_pow[i]);

        return reduce

    def reduce(self, xx: list[int]) -> list[int]:
        r = xx[:self.order]
        if len(r) < self.order:
            r += [0] * (self.order - len(r))
            return r
        for exp, v in enumerate(r[self.order:self.order*2], self.order):
            self.mac_in_place(r, v, self.small_pow[exp])
        power = self.small_pow[-1]
        for v in r[self.order*2:]:
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
        if self.order <= 1:
            return None

        me = list(self.c)
        me.append(1)
        derivative = me[1:]
        for i, x in enumerate(derivative):
            derivative[i] = x * (i+1) % self.p
        _, _, gcd = euclid(me, derivative, self.p)
        if len(gcd) > 1:
            return gcd

        # Try for gcd with X^{p^i} - X for 1 ≤ 2i ≤ order.
        X = self.small_pow[1];
        power = X                       # X¹ = X^p⁰
        for i in range(1, (self.order + 2) >> 1):
            power = self.pow(power, self.p)
            poly = list(power)
            # FIXME poly == 0 is possible here, if we split into polynomials
            # all of order i.
            poly[1] = (poly[1] - 1) % self.p
            _, _, gcd = euclid(me, poly, self.p)
            if len(gcd) != 1:
                if len(gcd) < self.order:
                    return gcd
                else:
                    assert self.order % i == 0
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

def poly_mult_scalar_in_place(xx: list[int], v: int, p: int):
    for i, x in enumerate(xx):
        xx[i] = x * v % p

def poly_div_mod(rr: list[int], yy: list[int], p: int) -> list[int]:
    # rr is reduced in place, quotient is returned.
    order = len(yy) - 1
    if order < 0:
        raise ZeroDivisionError('polynomial division by zero')
    while yy[order] == 0:
        if order == 0:
            raise ZeroDivisionError('polynomial division by zero')
        order -= 1
    if len(rr) <= order:
        return []
    monica = pow(yy[order], -1, p)      # Possible div-by-zero!
    result = [0] * (len(rr) - order)
    for i in reversed(range(len(result))):
        f = rr[i + order] * monica % p
        result[i] = f
        rr[i + order] = 0
        for j in range(order):
            rr[i + j] = (rr[i + j] - f * yy[j]) % p
    return result

def trim(xx: list[int]) -> list[int]:
    if xx == [] or xx[-1]:
        return xx
    for order in reversed(range(len(xx) - 1)):
        if xx[order]:
            return xx[:order+1]
    return []

def euclid(xx: list[int], yy: list[int], p: int) -> Tuple[
        list[int], list[int], list[int]]:
    aa = trim(list(xx))
    bb = trim(list(yy))
    k_a_x = [1]
    k_a_y: list[int] = []
    k_b_x: list[int] = []
    k_b_y = [1]
    while aa:
        prev_bb = list(bb)
        q = poly_div_mod(bb, aa, p)
        bb = trim(bb)
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

if __name__ == '__main__':
    a = [1,1,1,1,1,1]
    b = [1,1,1]
    assert poly_sub(a, b, 3) == [0,0,0,1,1,1]
    assert poly_sub(b, a, 3) == [0,0,0,2,2,2]

    print(poly_mult([1,1,1], [1,-1,1], 5))
    m = Monic(5, [1, 0, 1, 0]);
    assert m.mult([1,1,1],[1,-1,1]) == [0,0,0,0]
    assert m.crack() == 2               # Splits into 'unknown' order 2 factors.

    assert Monic(7, [4,3,4]).crack() is None
    assert Monic(7, [6,6,6,4]).crack() is None
    print(poly_mult([4,3,4,1], [6,6,6,4,1], 7))
    m = Monic(7, [3,0,3,1,4,4,1])
    assert m.mult([4,3,4,1,0,0,0], [6,6,6,4,1,0,0]) == [0] * 7
    assert m.crack() in ([4,3,4,1], [6,6,6,4,1])
