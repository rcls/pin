# An extremely naive quadratic sieve implementation.  The cubic gaussian
# elimination is awful!  The sieving has no enhancements.
import misc, quad

import bisect, math
from array import array
from typing import MutableSequence, Sequence
from dataclasses import dataclass

import math

class FoundFactor(Exception): pass

@dataclass
class Relation:
    sqrt: int
    exponents: Sequence[int]

@dataclass(slots=True)
class Sieve:
    N: int
    base: Sequence[int]
    log16: int
    relations: list[Relation]
    B: int
    bsize: int
    small_factors: list[int]
    factor_work: int
    blocks: list['SieveBlock']

    def __init__(self, N: int, maxB: int = 1<<30, k=None) -> None:
        self.N = N
        sqrtN, _ = misc.floor_sqrt(N)
        sqrtN = sqrtN + 1               # For convenience.
        self.relations = []
        twos = max(0, N.bit_length() - 32)
        self.log16 = int(round(math.log2(N) * 16))
        logN = math.log(N)
        B = int(math.exp(math.sqrt(logN * math.log(logN) / 2)))
        B = max(100, B)
        self.factor_work = max(100, int(math.sqrt(B) * 50))
        self.bsize = min(max(B, 1024), 1048576)
        print(self.bsize, min(maxB, B), B, sqrtN)
        B = min(maxB, B)
        self.B = B
        self.small_factors = [2, 0] * ((B + 1) // 2)
        for p in misc.sieve_primes_to(B):
            self.small_factors[p] = p
            for i in range(p * p, B, p):
                self.small_factors[i] = p
        self.blocks = []

        if k is None:
            #k = max((1,2,3,5,6,7,10,11,13,14,15,17,19), key=self.quality)
            k = max((1,2,3,5,6,7,10), key=self.quality)
        self.blocks.append(SieveBlock(self, k))

        pps = set((-1,2,))
        for b in self.blocks:
            pps.update(self.small_factors[pp] for pp in b.strides)
        self.base = list(pps)
        self.base.sort()
        print(f'{k=}, base len:', len(self.base))

    def quality(self, k: int) -> float:
        kN = self.N * k
        if k % 2 == 0:
            q = 0.5
        else:
            q = 1.0
        for p in misc.small_primes:
            if p == 2:
                continue
            j = misc.jacobi(kN, p)
            if j == 1:
                q += 2 / (p - 1)
            elif j == 0:
                q += 1 / p
        print(f'{k=} {q=}')
        return q

    def sieve(self) -> None:
        while len(self.relations) < len(self.base) + 10:
            sb = max(self.blocks, key=lambda b: b.last_num_rels)
            sb.sieve_block()
            print('Relations', len(s.relations), 'base len', len(s.base))

    def add_relation(self, x: int, block: 'SieveBlock') -> None:
        smooth = x * x - block.kN
        exps = [0 for _ in self.base]
        if smooth < 0:
            exps[0] = 1
        self.add_factors(abs(smooth), exps)
        self.relations.append(Relation(x, exps))

        # Back check.
        product = 1
        for i, p in enumerate(self.base):
            product *= pow(p, exps[i])
        assert smooth == product

    def add_factors(self, remain: int, exps: list[int]) -> None:
        while remain > self.B:
            factor = self.pollard_rho(remain)
            assert factor is not None, f'{remain}'
            remain //= factor
            if remain < factor:
                remain, factor = factor, remain
            self.add_factors(factor, exps)
        while remain > 1:
            p = self.small_factors[remain]
            assert p, f'{remain}'
            remain //= p
            index = bisect.bisect_left(self.base, p)
            assert 0 <= index < len(self.base), f'{remain} {p} {index}'
            assert self.base[index] == p, f'{remain} {p}'
            exps[index] += 1

    def pollard_rho(self, n: int) -> int|None:
        #print(f'pollard_rho', n)
        count = self.factor_work
        for i in range(5, n, 2):
            if n % i == 0:
                return i
            slow = i
            fast = i
            while True:
                slow = (slow * slow) % n + 1
                fast = (fast * fast) % n + 1
                fast = (fast * fast) % n + 1
                if slow == fast:
                    break
                g = math.gcd(slow - fast, n)
                if g != 1:
                    return g
                count -= 1
                if count == 0:
                    return None
        assert False

    def eliminate(self) -> None:
        # Build an array of byte arrays.
        rels = len(self.relations)
        bases = len(self.base)
        matrix = [0] * bases
        for i, rel in enumerate(self.relations):
            bit = 1 << i
            for j, e in enumerate(rel.exponents):
                if e & 1:
                    matrix[j] |= bit
        # Which rows have been assigned to select which relation.
        assigned: list[None|int] = [None] * bases
        for col in range(rels):
            # Try and deal with this column.
            bit = 1 << col
            non_zero = False
            for row, the_row in enumerate(matrix):
                if the_row & bit:
                    non_zero = True
                    if assigned[row] is None:
                        break
            else:                       # Failed to find this column.
                # Every set bit in this column is assigned.  Try it...
                if non_zero:
                    self.try_elimination(matrix, assigned, col, bit)
                continue
            assigned[row] = col
            for r, the_r in enumerate(matrix):
                if r != row and the_r & bit:
                    matrix[r] ^= the_row

        assert False, 'Damn'

    def try_elimination(self, matrix: list[int],
                        assigned: list[int|None], col: int, bit: int) -> None:
        column_to_null = [1 if row & bit else 0 for row in matrix]
        columns = [col]
        for row, e in enumerate(column_to_null):
            if e:
                c = assigned[row]
                assert c is not None
                columns.append(c)
        # Now construct the x and the exponents.
        exps = [0] * len(self.base)
        x = 1
        for c in columns:
            rel = self.relations[c]
            x = x * rel.sqrt % self.N
            for i, e in enumerate(rel.exponents):
                exps[i] += e
        y = 1
        for i, e in enumerate(exps):
            assert e & 1 == 0
            y = y * pow(self.base[i], e // 2) % self.N
        assert x * x % self.N == y * y % self.N
        f = math.gcd(self.N, x - y)
        if f == self.N:
            f = math.gcd(self.N, x + y)
        if 1 < f < self.N:
            raise FoundFactor(self.N, f)

class SieveBlock:
    S: Sieve
    k: int
    kN: int
    pos_offset: MutableSequence[int]
    neg_offset: MutableSequence[int]
    strides: Sequence[int]
    logs: Sequence[int]
    pos_start: int
    neg_start: int
    last_num_rels: int
    def __init__(self, S: Sieve, k: int):
        self.S = S
        self.k = k
        self.kN = k * S.N
        sqrtkN, _ = misc.floor_sqrt(self.kN)
        sqrtkN += 1
        self.pos_start = sqrtkN
        self.neg_start = sqrtkN

        self.strides = array('I')
        self.logs = array('I')
        self.pos_offset = array('I')
        self.neg_offset = array('I')
        self.last_num_rels = S.bsize

        #quality = 0.0
        for p in misc.sieve_primes_to(S.B):
            if p == 2:
                continue
            #quality += misc.jacobi(self.kN, p) / p
            # p|k is a special case.
            if k % p == 0:
                assert k // p % p != 0  # We only deal with square free k.
                if S.N % p == 0:
                    raise FoundFactor(S.N, p)
                sqrtkN_mod_p = sqrtkN % p
                pos = -sqrtkN_mod_p % p
                self.strides.append(p)
                self.pos_offset.append(pos)
                self.neg_offset.append(p - pos)
                self.logs.append(int(round(16 * math.log2(p))))
                continue

            ring_sqrt = quad.QuadRing(p).maybe_sqrt(self.kN % p)
            if ring_sqrt is None:
                continue
            if ring_sqrt == 0:
                raise FoundFactor(S.N, p)

            inv_2sqrt = pow(2 * ring_sqrt, -1, p)
            prev_log16pp = 0 #int(round(16 * math.log2(p)))
            pp = p
            while pp < S.B:
                # Append twice, once for each sqrt.
                self.strides.append(pp)
                self.strides.append(pp)

                log16pp = int(round(16 * math.log2(pp)))

                # Note that we only record log(p) as we process the lessor prime
                # powers, and are only an increment of one over that.
                self.logs.append(log16pp - prev_log16pp)
                self.logs.append(log16pp - prev_log16pp)
                prev_log16pp = log16pp

                sqrtkN_mod_pp = sqrtkN % pp
                # We will be looking for (pos_start + x)² ≡ N (mod pp).  So we
                # want pos s.t. sqrtN_mod_pp + pos ≡ ±ring_sqrt
                pos = (ring_sqrt - sqrtkN_mod_pp) % pp
                self.pos_offset.append(pos)
                self.neg_offset.append(pp - pos) # Note that we skip zero.
                # Repeat but change the sign of ring_sqrt.
                pos = (-ring_sqrt - sqrtkN_mod_pp) % pp
                self.pos_offset.append(pos)
                self.neg_offset.append(pp - pos) # Note that we skip zero.

                next_pp = pp * p
                # Lift ring_sqrt to mod next_pp
                # We want ring_sqrt² ≡ N (mod next_pp).
                kN_mod_next_pp = self.kN % next_pp
                adjust = (kN_mod_next_pp - ring_sqrt * ring_sqrt) % next_pp
                assert adjust % pp == 0
                # ring_sqrt += adjust / f'(ring_sqrt)
                # where f' is the derivative of x² mod p, i.e., 2*ring_sqrt.
                ring_sqrt = ring_sqrt + adjust * inv_2sqrt
                assert ring_sqrt * ring_sqrt % next_pp == kN_mod_next_pp
                pp = next_pp

    def sieve_block(self) -> None:
        S = self.S
        pos = array('I')
        neg = array('I')
        # Deal to powers of two.
        for i in range(S.bsize):
            x = (self.pos_start & 0xffffffff) + i
            xxkN = x * x - (self.kN & 0xffffffff)
            xlow = xxkN & -xxkN
            if xlow == 0:
                pos.append(32 * 16)
            else:
                pos.append((xlow.bit_length() - 1) * 16)

            y = (self.neg_start & 0xffffffff) - i
            yykN = y * y - (self.kN & 0xffffffff)
            ylow = yykN & -yykN
            if ylow == 0:
                neg.append(32 * 16)
            else:
                neg.append((ylow.bit_length() - 1) * 16)

        for i, p in enumerate(self.strides):
            l = self.logs[i]
            j = self.pos_offset[i]
            while j < S.bsize:
                pos[j] += l
                j += p
            self.pos_offset[i] = j - S.bsize

            j = self.neg_offset[i]
            while j < S.bsize:
                neg[j] += l
                j += p
            self.neg_offset[i] = j - S.bsize

        pos_xxn = self.pos_start * self.pos_start - self.kN
        neg_xxn = self.neg_start * self.neg_start - self.kN
        num_rels = len(S.relations)
        for i in range(S.bsize):
            pos_error = math.log2(pos_xxn) * 16 - pos[i]
            if -16 < pos_error < 16:
                S.add_relation(self.pos_start + i, self)

            neg_error = math.log2(abs(neg_xxn)) * 16 - neg[i]
            if -16 < neg_error < 16:
                S.add_relation(self.neg_start - i, self)

            # x = self.pos_start + i
            pos_xxn += 1 + 2 * (self.pos_start + i)
            neg_xxn += 1 - 2 * (self.neg_start - i)
        self.last_num_rels = len(S.relations) - num_rels

        self.pos_start += S.bsize
        self.neg_start -= S.bsize

if __name__ == '__main__':
    import time
    import sys
    #s = Sieve(549772066901 * 1099524735011)

    start = time.time()

    # 42 secs, B=100000
    #s = Sieve(9223372088242534787 * 18446744556051857693, maxB = 100000)

    # 711 secs, B=200000 [1].
    # 150 bits.
    # 500000, base len 20836, sieve 391, elim 231, total 622
    # 400000, base len 16923, sieve 409, elim 145, total 556
    # 300000, base len 13028, sieve 452, elim 78, total 530
    # 500000k2, base len 20803, sieve 250, elim 209, total 459
    # 400000k2, base len 16966, sieve 235, elim 133, total 369
    # 300000k2, base len 13055, sieve 245, elim 72, total 317
    # 250000k2, base len 11087, sieve 256, elim 53, total 309
    # 230000k2, base len 10295, sieve 260, elim 40, total 300
    # 220000k2, base len 9870, sieve 269, elim 37, total 306
    # 200000k2, base len 9052, sieve 280, elim 34, total 314

    #*250000k5, base len 11087, sieve 178, elim 48, total 226
    # 230000k5, base len 10307, sieve 184, elim 44, total 229
    #s = Sieve(18889465942126000625021 * 37778931879853556277577, maxB=250000)

    # 156 bits
    # 1000000, base len 39107, sieving 408, elim 1051, total 1460.
    #  500000, base len 20640, sieve 224, elim 212, total 438.
    #  300000, base len 12878, sieve 218, elim 73, total 291.
    #  290000, base len 12490, sieve 216, elim 67, total 284
    #  280000, base len 12097, sieve 221, elim 58, total 280
    #  270000, base len 11694, sieve 219, elim 57, total 276
    #  260000, base len 11298, sieve 222, elim 50, total 272
    #  250000, base len 10929, sieve 228, elim 47, total 275
    #  240000, base len 10531, sieve 231, elim 43, total 275
    #  230000, base len 10139, sieve 232, elim 43, total 275
    #s = Sieve(151115727451968613138879 * 302231454942896865282311, maxB=int(sys.argv[1]), k = 1)

    # 180 bits, k=1
    # 600000, base len 24473
    # 1000000, base len 39159, sieve 2332, elim 1245, total 3579
    #  900000, base len 35581, sieve 2402, elim  983, total 3386
    #  800000, base len 31905, sieve 2539, elim  729, total 3270
    #  700000, base len 28165, sieve 2660. elim  551, total 3213
    #  600000, 24473, 2875, 384, 3260
    s = Sieve(618970019644176682696314337 * 1237940039291550758926781041, k=1,
              maxB = int(sys.argv[1]))

    for sb in s.blocks:
        assert len(sb.strides) == len(sb.pos_offset)
        assert len(sb.strides) == len(sb.neg_offset)
        assert len(sb.strides) == len(sb.logs)
        for i in range(2, len(sb.strides)):
            p = sb.strides[i]
            pos = sb.pos_start + sb.pos_offset[i]
            assert (pos * pos - sb.kN) % p == 0, f'{p}'
            neg = sb.neg_start - sb.neg_offset[i]
            assert (neg * neg - sb.kN) % p == 0

    start_sieve = time.time()
    print('Setup took', start_sieve - start)
    print('Sieving...')
    s.sieve()
    sieve_end = time.time()
    print('Sieving took', sieve_end - start_sieve)
    try:
        s.eliminate()
    except FoundFactor as e:
        print('Factor', e.args[1])
    end = time.time()
    print('Elimination took', end - sieve_end)
    print('Total', end - start, 'B', s.B)
