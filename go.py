#!/usr/bin/python3

from misc import small_primes
from pollard_rho import pollard_rho, unique_prime_factors
from pratt_cert import pratt_cert

import argparse

parse = argparse.ArgumentParser()
parse.add_argument('--certify', '-c', action='store_true')
parse.add_argument('--partial', '-p', action='store_true')
parse.add_argument('expr', nargs='*')
args = parse.parse_args()

def process(v, s):
    if not args.certify and not args.partial:
        print(f'{s}:', ' '.join(str(x) for x in unique_prime_factors(v)))
        return
    c = pratt_cert(v, {})
    if c is not None:
        c.verify()
        print(f'{s} is prime')
    elif v > 1 and args.partial:
        f = pollard_rho(v)
        print(f'{s} has factor {f}')
    else:
        print(f'{s} is NOT prime')

def is_moderate_prime(n):
    assert n < 400 * 400
    if n in small_primes or n < 2:
        return True
    for p in small_primes:
        if n % p == 0:
            return False
        return True

def F(n: int) -> int:
    return (1 << (1 << n)) + 1
def M(n: int) -> int:
    return (1 << n) - 1

for s in args.expr:
    process(eval(s), s)

if len(args.expr) < 1:
    for s in range(400 * 400):
        if is_moderate_prime(s):
            print(f'--------- M({s}) ---------')
            process(M(s), f'M({s})')

