
from math import gcd, lcm
from typing import Tuple
from numbers import Rational
from fractions import Fraction

def rational(f: Rational) -> str:
    n, d = f.numerator, f.denominator

    if d == 1:
        return str(n)

    s = '-' if n < 0 else ''
    an = abs(n)
    if an == 1 and d <= 10:
        return s + '01½⅓¼⅕⅙⅐⅛⅑⅒'[d]
    if an + 1 == d and d < 7:
        return '0½⅔¾⅘⅚'[an]
    if d == 5 and an < 5:
        return s + '0⅕⅖⅗⅘'[an]
    if d == 8 and an < 8:
        return s + '0⅛¼⅜½⅝¾⅞'[an]
    return f'{n}/{d}'

def ratlin_basic(a: Rational, b: Rational, k: str) -> str:
    if b == 0:
        return rational(a)

    if abs(b.numerator) != 1:
        bb = rational(Fraction(abs(b.numerator), abs(b.denominator)))
        if '/' in bb:
            bbk = bb + ' × ' + k
        else:
            bbk = bb + ' ' + k
    elif b.denominator == 1:
        bbk = k
    else:
        bbk = k + '/' + str(b.denominator)

    if a == 0:
        return '-' + bbk if b < 0 else bbk

    if a.numerator < 0 and b.numerator >= 0:
        return bbk + ' - ' + rational(-a)

    return rational(a) + (' + ' if b.numerator >= 0 else ' - ') + bbk

def factorout(o: Rational, a: Rational, b: Rational, k: str) -> str:
    inner = ratlin_basic(a / o, b / o, k)
    if '+' in inner or '-' in inner:
        inner = '(' + inner + ')'

    if o.numerator == 1:
        return inner + '/' + str(o.denominator)
    if o.numerator == -1:
        return '-' + inner + '/' + str(o.denominator)

    return rational(o) + ' ' + inner

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
