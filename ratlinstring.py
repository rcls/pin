
from math import gcd, lcm
from typing import Tuple
from numbers import Rational
from fractions import Fraction

def simple(n: int, d: int) -> str|None:
    if n == 1 and d <= 10:
        return '01½⅓¼⅕⅙⅐⅛⅑⅒'[d]
    if n + 1 == d and 2 <= d <= 6:
        return '0½⅔¾⅘⅚'[n]
    if d == 5 and n < 5:
        return '0⅕⅖⅗⅘'[n]
    if d == 8 and n < 8:
        return '0⅛¼⅜½⅝¾⅞'[n]
    if d == 1:
        return str(n)
    return None

def rational(f: Rational, k:str = '') -> str:
    n, d = f.numerator, f.denominator
    sign = '-' if n < 0 else ''
    an = abs(n)

    if an == 0:
        return '0'
    if an == 1 and d == 1 and k != '':
        return sign + k.strip()
    if d == 1:
        return str(n) + k.rstrip()

    s = simple(an, d)
    if s is not None:
        return sign + s + k.rstrip()
    if k.endswith(' '):
        sp = ' '
    else:
        sp = ''
    k = k.strip()
    if an == 1 and k != '':
        return f'{sign}{k}{sp}/{sp}{d}'.lstrip()
    return f'{n}{sp}{k}{sp}/{sp}{d}'

def ratlin_basic(a: Rational, b: Rational, k: str) -> str:
    if b == 0:
        return rational(a)

    bk = rational(Fraction(abs(b.numerator), abs(b.denominator)), k)
    if a == 0:
        return '-' + bk if b < 0 else bk

    if a.numerator < 0 and b.numerator >= 0:
        return bk + ' - ' + rational(-a)

    return rational(a) + (' + ' if b.numerator >= 0 else ' - ') + bk

def factorout(o: Rational, a: Rational, b: Rational, k: str) -> str:
    inner = ratlin_basic(a / o, b / o, k)
    if '+' in inner or '-' in inner:
        inner = '(' + inner + ') '
    else:
        inner = ' ' + inner + ' '
    return rational(o, inner)

def ratlinstring(a: Rational, b: Rational, k: str) -> str:
    # Format a + b × k in the most economical way.
    # Try taking out various factors:
    an, ad = abs(a.numerator), a.denominator
    bn, bd = abs(b.numerator), b.denominator
    result = ratlin_basic(a, b, k)
    for on in lcm(an, bn), max(an, bn), min(an, bn), gcd(an, bn), 1:
        for od in lcm(ad, bd), max(ad, bd), min(ad, bd), gcd(ad, bd), 1:
            if on != 0 and od != 0:
                for son in on, -on:
                    attempt = factorout(Fraction(son, od), a, b, k)
                    if len(attempt) < len(result):
                        result = attempt
    return result

if __name__ == '__main__':
    from sys import argv
    if len(argv) == 2:
        print(ratlinstring(Fraction(argv[1]), Fraction(argv[1]), 'φ'))
    elif len(argv) == 3:
        print(ratlinstring(Fraction(argv[1]), Fraction(argv[2]), 'φ'))
    elif len(argv) == 3:
        print(ratlinstring(Fraction(argv[1]), Fraction(argv[2]), argv[3]))
    else:
        print('One to three arguments please')
