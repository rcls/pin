import math, time

delta_by_bits = {
    32: 19, 33: 40, 34: 11, 35: 15, 36: 16, 37: 108, 38: 45, 39: 24,
    40: 16, 41: 7, 42: 15, 43: 12, 44: 18,
    45: 9, 46: 9, 47: 7, 48: 2, 49: 59,
    50: 60, 51: 76, 52: 4, 53: 10, 54: 17,
    55: 19, 56: 24, 57: 67, 58: 20, 59: 9,
    60: 34, 61: 16, 62: 12, 63: 43, 64: 12,
    65: 57, 66: 37, 67: 26, 68: 99, 69: 19,
    70: 28, 71: 50, 72: 30, 73: 59, 74: 10,
    75: 39, 76: 62, 77: 192, 78: 1, 79: 36,
    80: 212, 81: 1, 82: 22, 83: 23, 84: 7,
    85: 37, 86: 44, 87: 44, 88: 76, 89: 94,
    90: 43, 91: 88, 92: 26, 93: 21, 94: 81,
    95: 23, 96: 71, 97: 106, 98: 192, 99: 133,
    100: 49, 101: 88,

    128: 55, 129: 193, 130: 66, 256: 253, 257: 154, 258: 4,
    512: 130, 513: 369, 514: 264, 1024: 903, 1025: 301, 1026: 666,
    2048: 1280, 2049: 547, 2050: 3845, 4096: 1915, 4097: 977, 4098: 8563,
    8192: 1625, 8193: 14319, 8194: 52, 16384: 25567, 16385: 1086, 16386: 17472,
    32768: 9274, 32769: 9033, 32770: 14620,
    65536: 59525, 65537: 108747}

prime_by_bits = {
    2: 3, 3: 5, 4: 11, 5: 17, 6: 37, 7: 67, 8: 131, 9: 257, 10: 521,
    11: 1031, 12: 2053, 13: 4099, 14: 8209, 15: 16411, 16: 32771, 17: 65537,
    18: 131101, 19: 262147, 20: 524309, 21: 1048583, 22: 2097169,
    23: 4194319, 24: 8388617, 25: 16777259, 26: 33554467, 27: 67108879,
    28: 134217757, 29: 268435459, 30: 536870923, 31: 1073741827,
}

def _populate() -> None:
    for bits, delta in delta_by_bits.items():
        subbits = (bits + 3) // 2
        p = prime_by_bits[subbits]
        start_t = (1 << bits-1) // p
        t = start_t + delta
        prime = p * t + 1
        assert prime & 1 == 1, f'{bits=}'
        assert prime.bit_length() == bits
        prime_by_bits[bits] = prime

_populate()

prime4096 = prime_by_bits[4096]
prime8192 = prime_by_bits[8192]
prime16384 = prime_by_bits[16384]
prime32768 = prime_by_bits[32768]
prime65536 = prime_by_bits[65536]

def test_len() -> None:
    assert prime4096.bit_length() == 4096
    assert prime8192.bit_length() == 8192
    assert prime16384.bit_length() == 16384
    assert prime32768.bit_length() == 32768
    assert prime65536.bit_length() == 65536
    for bits, prime in prime_by_bits.items():
        assert prime.bit_length() == bits

def verify_one_fermat(bits: int, base: int, verbose:bool = False) -> None:
    start = time.time()
    prime = prime_by_bits[bits]
    if verbose:
        print(f'Test Fermat {bits}')
    if bits in delta_by_bits:
        assert pow(base, prime-1, prime) == 1
    else:
        # Use a Pratt cert instead.
        import pratt_cert
        assert bits < 100
        assert pratt_cert.pratt_cert(prime, {})
    if verbose:
        print(f'Passed Fermat {bits} in', time.time() - start, 'seconds')

def verify_one_sub(bits: int, base: int, verbose:bool = False) -> None:
    if not bits in delta_by_bits:
        return
    if verbose:
        print(f'Test sub {bits}')
        start = time.time()
    prime = prime_by_bits[bits]
    subbits = (bits + 3) // 2
    p = prime_by_bits[subbits]
    start_t = (1 << bits-1) // p
    t = start_t + delta_by_bits[bits]
    assert t < p
    assert prime == p * t + 1
    assert math.gcd(pow(base, t, prime) - 1, prime) == 1
    if verbose:
        print(f'Passed sub {bits} in', time.time() - start, 'seconds')

def test_verify_serial() -> None:
    import pseudo_prime
    for bits in prime_by_bits:
        if bits < 8000:
            verify_one_fermat(bits, 3)
            verify_one_sub(bits, 3)
            assert pseudo_prime.baillie_psw(prime_by_bits[bits])

def test_verify_parallel(limit: int = 30000) -> None:
    import joblib
    joblib.Parallel(n_jobs=-1, batch_size=1)(
        joblib.delayed(func)(bits, 3, True)
        for bits in reversed(prime_by_bits) if 8000 <= bits < limit
        for func in (verify_one_fermat, verify_one_sub))

def regenerate(bits: int, low16: int, wraps: int = 0) -> None:
    subbits = (bits + 3) // 2
    p = prime_by_bits[subbits]
    start_t = (1 << bits-1) // p
    t = (start_t & -65536) + low16
    if t < start_t:
        t += 65536
    t += wraps * 65536
    print('Delta', bits, t - start_t)
    prime = p * t + 1
    assert pow(2, prime - 1, prime) == 1
    assert math.gcd(pow(2, t, prime) - 1, prime) == 1

def generate(bits: int) -> None:
    #if bits in prime_by_bits:
    #    return
    import pocklington
    prime = pocklington.pocklington_serial(bits)
    subbits = (bits + 3) // 2
    p = prime_by_bits[subbits]
    start_t = (1 << bits-1) // p
    assert prime % p == 1
    t = prime // p
    assert t < p
    assert pow(2, prime-1, prime) == 1
    assert math.gcd(pow(2, t, prime) - 1, prime) == 1
    print(f'Delta for {bits} is {t - start_t}')
    delta_by_bits[bits] = t - start_t
    prime_by_bits[bits] = prime

if __name__ == '__main__':
    import sys
    regenerate(int(sys.argv[1]), int(sys.argv[2]))
