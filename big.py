import sys

delta_dict = {
    128: 55, 129: 193, 130: 66, 256: 253, 257: 154, 258: 4, 512: 130, 513: 369,
    514: 264, 1024: 903, 1025: 301, 1026: 666, 2048: 1280, 2049: 547,
    2050: 3845, 4096: 1915, 4097: 977, 4098: 8563, 8192: 1625, 8193: 14319,
    8194: 52, 16384: 25567, 16385: 1086, 16386: 17472, 32768: 9274, 32769: 9033,
    65536: 59525}

prime_by_bits = {
    32: 2148696083, 64: 9223372088242534787, 65: 18446744556051857693,
    66: 36893488459164099677}

def _populate():
    for bits, delta in delta_dict.items():
        subbits = (bits + 3) // 2
        p = prime_by_bits[subbits]
        start_t = (1 << bits-1) // p
        t = start_t + delta
        assert t & 1 == 0
        prime = p * t + 1
        assert prime.bit_length() == bits
        prime_by_bits[bits] = prime

_populate()

prime4096 = prime_by_bits[4096]
prime8192 = prime_by_bits[8192]
prime16384 = prime_by_bits[16384]
prime32768 = prime_by_bits[32768]
prime65536 = prime_by_bits[65536]

bigs = prime4096, prime8192, prime16384, prime32768, prime65536

def test_len():
    assert prime4096.bit_length() == 4096
    assert prime8192.bit_length() == 8192
    assert prime16384.bit_length() == 16384
    assert prime32768.bit_length() == 32768
    assert prime65536.bit_length() == 65536

def test4096():
    import pseudo_prime
    assert pseudo_prime.baillie_psw(prime4096)

def test8192():
    import pseudo_prime
    assert pseudo_prime.baillie_psw(prime8192)

#def test16384():
#    import pseudo_prime
#    assert pseudo_prime.miller_rabin1(prime16384, 3)
    #assert pseudo_prime.baillie_psw(prime16384)

#def test32768():
#    # Hard coded M-R test.
#    assert prime32768 & 7 == 5
#    assert pow(-3, (prime32768-1) // 4, prime32768) == 1

def verify_one_fermat(bits, base):
    prime = prime_by_bits[bits]
    if not bits in delta_dict:
        assert bits < 100
        print(f'Pratt cert check {bits}')
        import pratt_cert
        assert pratt_cert.pratt_cert(prime, {})
        return
    import time
    start = time.time()
    print(f'Fermat {bits}')
    assert pow(base, prime-1, prime) == 1
    print(f'Passed fermat {bits} in', time.time() - start, 'seconds')

def verify_one_sub(bits, base):
    from math import gcd
    if not bits in delta_dict:
        return
    print(f'Sub test {bits}')
    import time
    start = time.time()
    prime = prime_by_bits[bits]
    subbits = (bits + 3) // 2
    p = prime_by_bits[subbits]
    start_t = (1 << bits-1) // p
    t = start_t + delta_dict[bits]
    assert prime == p * t + 1
    assert gcd(pow(base, t, prime) - 1, prime) == 1
    print(f'Passed sub test {bits} in', time.time() - start, 'seconds')

def verify_one(bits, base, which=None):
    if which is None or which == 1:
        verify_one_fermat(bits, base)
    if which is None or which == 2:
        verify_one_sub(bits, base)

def verify(base=2):
    for bits in prime_by_bits:
        verify_one(bits, base)

def test_verify_small():
    from joblib import Parallel, delayed
    for bits in prime_by_bits:
        if bits < 8000:
            verify_one(bits, 2)
    Parallel(n_jobs=-1, batch_size=1)(
        delayed(verify_one)(bits, 2, which)
        for which in (1, 2)
        for bits in reversed(prime_by_bits) if 8000 < bits < 30000)

def regenerate(bits, low16):
    subbits = (bits + 3) // 2
    p = prime_by_bits[subbits]
    start_t = (1 << bits-1) // p
    t = (start_t & -65536) | low16
    if t < start_t:
        t += 65536
    print('Delta', bits, t - start_t)
    prime = p * t + 1
    assert pow(2, prime - 1, prime) == 1
    import math
    assert math.gcd(pow(2, t, prime) - 1, prime) == 1
