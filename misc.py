
from math import sqrt

tiny_primes = 2, 3, 5, 7, 11, 13, 17, 19
small_primes = {
    x for x in range(2, 400)
    if x in tiny_primes or not any(x % p == 0 for p in tiny_primes)}
small_primes_list = sorted(small_primes)
