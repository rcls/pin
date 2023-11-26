
from dataclasses import dataclass

@dataclass
class Mersenne:
    p: int
    N: int
    def __init__(self, p: int):
        self.p = p
        self.N = (1 << p) - 1

    def mod(self, x: int) -> int:
        # Not general purpose, just enough to deal with multiplication.
        xx = abs(x)
        if xx > self.N:
            xx = (xx & self.N) + (xx >> self.p)
        if x < 0:
            xx = self.N - xx
        if xx.bit_length() >= self.p:
            return xx - self.N
        else:
            return xx

    def power(self, b: int, m: int) -> int:
        if m == 0:
            return 1
        result = b
        for i in reversed(range(m.bit_length() - 1)):
            result = self.mod(result * result)
            if m & (1 << i):
                result = self.mod(result * b)
        return result

    def S(self, m: int) -> int:
        s = 4
        for _ in range(m):
            s = self.mod(s * s) - 2
        return s

    def is_prime(self) -> bool:
        # Let n = 2ᵖ-1.  Work in ℤₙ[√3].  Jacobi (3|n) = -(3|n) = -1 as p is
        # odd.  So if n is prime, then 3 is not a q.r. & the extension is a
        # quadratic field.
        #
        # 2 has a square root mod n, since 2ᵖ⁺¹ ≡ 2 and p+1 is even.
        #
        # Let 𝛼 = (1 + √3)/√2, so that 𝛼◌̅𝛼 = -1.  Let 𝜔 = 𝛼² = 2 + √3, so that
        # 𝜔◌̅𝜔 = 1 and 𝜔+◌̅𝜔 = 4.
        #
        # By induction, S(k) ≡ 𝜔^(2ᵏ) + ◌̅𝜔^(2ᵏ) = 𝛼^(2ᵏ⁺¹) + ◌̅𝛼^(2ᵏ⁺¹) for all
        # k.
        #
        # If n is prime, then using Frobenius, 𝛼^(2ᵖ) = 𝛼ⁿ·𝛼 = ◌̅𝛼·𝛼 = -1, and
        # so:
        #
        # S(p-2) = 𝛼^(2ᵖ⁻¹) + ◌̅𝛼^(2ᵖ⁻¹) = ◌̅𝛼^(2ᵖ⁻¹)·(𝛼^(2ᵖ) + 1) = 0
        #
        # Conversely, suppose that S(p-2) ≡ 0 (mod n) but n is composite.
        #
        # Let q be the smallest prime factor of n, so that q² < 2ᵖ.
        #
        # Now, S(p-2) ≡ 0 (mod q) also.  Work in the ring ℤ_q[√3], and let 𝜔 =
        # 2 + √3 as before.
        #
        # 𝜔^(2ᵖ⁻²) + ◌̅𝜔^(2ᵖ⁻²) ≡ S(p-2) ≡ 0, so 𝜔^(2ᵖ⁻²) = -◌̅𝜔^(2ᵖ⁻²).
        #
        # Using 𝜔 = 1/◌̅𝜔 we get
        #
        # 𝜔^(2ᵖ⁻¹) = 𝜔^(2ᵖ⁻²) / ◌̅𝜔^(2ᵖ⁻²) = -1,
        #
        # so that the multiplicative mod q order of 𝜔 is 2ᵖ > q² ≥ #ℤ_q[√3],
        # which is impossible.
        return self.p == 2 or self.S(self.p - 2) == 0

mersennes = (
    2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281,
    3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497,
    86243, 110503, 132049, 216091, 756839, 859433, 1257787, 1398269, 2976221,
    3021377, 6972593, 13466917, 20996011, 24036583, 25964951, 30402457,
    32582657, 37156667, 42643801, 43112609, 57885161)

def test_mersennes() -> None:
    import joblib
    def mprime(p: int) -> None:
        assert Mersenne(p).is_prime()
    not_too_big = [p for p in mersennes if p < 40000]
    joblib.Parallel(n_jobs=-1, batch_size=1)(
        joblib.delayed(mprime)(p) for p in reversed(not_too_big))

def test_non_mersennes() -> None:
    import misc
    for p in misc.modest_primes_list[:250]:
        if not p in mersennes:
            assert not Mersenne(p).is_prime()

#import misc
#m11213 = Mersenne(11213)
#print(misc.timecall(m11213.is_prime))
#print(misc.timecall(m11213.power, 2, m11213.N-1))
#print(misc.timecall(pow, 2, m11213.N-1, m11213.N))
#
#print(misc.timecall(Mersenne(44497).is_prime))
