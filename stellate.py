from normalize import Quadratic
from fractions import Fraction
from dataclasses import dataclass
#from typing import TypeVar, Generic
from ratlinstring import ratlinstring

#T = TypeVar('T')

@dataclass
class Vector:
    x: Quadratic
    y: Quadratic
    z: Quadratic
    def __mul__(self, v: 'Vector') -> Quadratic:
        return self.x * v.x + self.y * v.y + self.z * v.z
    def __add__(self, v: 'Vector') -> 'Vector':
        return Vector(self.x + v.x, self.y + v.y, self.z + v.z)
    def scale(self, q: Quadratic) -> 'Vector':
        return Vector(self.x * q, self.y * q, self.z * q)
    def normsq(self) -> Quadratic:
        return self * self

def eliminate(A: list[list[Quadratic]], r: list[Quadratic]) -> None:
    # Solve A·x = r in-place.  Naïve elimination.
    size = len(A)
    assert len(r) == size
    assert min(len(v) for v in A) == size
    assert max(len(v) for v in A) == size
    base = A[0][0].b
    for v in A:
        for x in v:
            assert x.b == base
    for x in r:
        assert x.b == base
    for i in range(size):
        # Swap to being a non-zero element, if necessary:
        if not A[i][i]:
            for j in range(i + 1, size):
                if A[j][i]:
                    break
            else:
                #print(i, A[i][i], A, r)
                raise ZeroDivisionError('Singular matrix')
            A[i], A[j] = A[j], A[i]
            r[i], r[j] = r[j], r[i]
        inv = A[i][i]
        for j in range(i, size):
            A[i][j] /= inv
        r[i] /= inv
        for ii in range(size):
            if ii == i:
                continue
            m = A[ii][i]
            for j in range(i, size):
                A[ii][j] -= m * A[i][j]
            r[ii] -= m * r[i]

def lift(x: int | Fraction | Quadratic) -> Quadratic:
    if type(x) == Quadratic:
        return x
    if type(x) == Fraction:
        return Quadratic(x, Fraction(0), 5)
    assert type(x) == int
    return Quadratic(Fraction(x), Fraction(0), 5)

def liftv(v: list[int | Fraction | Quadratic]) -> list[Quadratic]:
    return [lift(x) for x in v]

def liftm(A: list[list[int | Fraction | Quadratic]]) -> list[list[Quadratic]]:
    return [liftv(v) for v in A]


def tri_intersect(a: list[Vector], b: list[Vector], c: list[Vector]) -> list[Quadratic]:
    r = liftv([1, 1, 1, 0, 0, 0, 0, 0, 0])
    A = liftm(
        [[1, 1, 1, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 1, 1, 1, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 1, 1, 1],
         [a[0].x, a[1].x, a[2].x, -b[0].x, -b[1].x, -b[2].x, 0, 0, 0],
         [a[0].y, a[1].y, a[2].y, -b[0].y, -b[1].y, -b[2].y, 0, 0, 0],
         [a[0].z, a[1].z, a[2].z, -b[0].z, -b[1].z, -b[2].z, 0, 0, 0],
         [0, 0, 0, b[0].x, b[1].x, b[2].x, -c[0].x, -c[1].x, -c[2].x],
         [0, 0, 0, b[0].y, b[1].y, b[2].y, -c[0].y, -c[1].y, -c[2].y],
         [0, 0, 0, b[0].z, b[1].z, b[2].z, -c[0].z, -c[1].z, -c[2].z],
         ])
    eliminate(A, r)
    return r

def rotv(a: Vector) -> Vector:
    return Vector(a.z, a.x, a.y)

gold = Quadratic(Fraction(1,2), Fraction(1,2), 5)
one = gold / gold
zero = one - one
ico_points = [Vector(zero, s, t) for s in (-one, one) for t in (-gold, gold)]

ico_points = ico_points + [rotv(x) for x in ico_points] \
    + [rotv(rotv(x)) for x in ico_points]

def rot5(v: Vector) -> Vector:
    half = Quadratic(Fraction(1,2), Fraction(0), 5)
    return Vector(
        half * v.x + half * gold * v.y + half / gold * v.z,
        -half * gold * v.x + half / gold * v.y + half * v.z,
        half/gold * v.x - half * v.y + half * gold * v.z)

def rot5_2(v: Vector) -> Vector:
    return rot5(rot5(v))
def rot5_3(v: Vector) -> Vector:
    return rot5_2(rot5(v))
def rot5_4(v: Vector) -> Vector:
    return rot5_2(rot5_2(v))
def rot5_5(v: Vector) -> Vector:
    return rot5(rot5_4(v))

vx = Vector(lift(1), lift(0), lift(0))
vy = rotv(vx)
vz = rotv(vy)
#print(vx, vy, rot5(vz))

for a in vx, vy, vz:
    assert rot5_5(a) == a
    assert rot5(a) != a
    for b in vx, vy, vz:
        assert rot5(a) * rot5(b) == (one if a == b else zero)

vertex = Vector(zero, one, gold)
tri = [vertex, rotv(vertex), rotv(rotv(vertex))]
def rx(p: Vector) -> Vector:
    return Vector(p.x, -p.y, -p.z)
def ry(p: Vector) -> Vector:
    return Vector(-p.x, p.y, -p.z)

def four(t: list[Vector]) -> list[list[Vector]]:
    return [t, [rx(p) for p in t], [ry(p) for p in t], [ry(rx(p)) for p in t]]

def five(t: list[Vector]) -> list[list[Vector]]:
    result = [t]
    for i in range(4):
        t = [rot5(x) for x in t]
        result.append(t)
    return result

def twenty(t: list[Vector]) -> list[list[Vector]]:
    return sum((five(u) for u in four(t)), [])

def goldrat(q: Quadratic) -> str:
    assert q.b == 5
    a = q.r - q.q
    b = 2 * q.q
    z = Fraction(0)
    assert q == Quadratic(a, z, 5) + Quadratic(b, z, 5) * gold
    return ratlinstring(a, b, 'φ')

if __name__ == '__main__':
    faces = twenty(tri)

    for i in range(1, len(faces)):
        for j in range(i+1, len(faces)):
            u, v = faces[i], faces[j]
            try:
                w = tri_intersect(tri, u, v)
                p = tri[0].scale(w[0]) + tri[1].scale(w[1]) + tri[2].scale(w[2])
                dd = p.normsq() / vertex.normsq()
                print(float(dd), goldrat(dd), dd, p)
            except ZeroDivisionError:
                pass
