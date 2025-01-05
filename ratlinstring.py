
from math import gcd, lcm
from typing import Tuple
from numbers import Rational
from fractions import Fraction

def simple(n: int, d: int) -> str|None:
    if n == 1 and d == 10:
        return '⅒ '
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

#⅕(175 + 78√5)
#35 + 78/5 √5
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
    if sp == '':
        return f'{n}{k}/{d}' #{sp}{k}{sp}/{sp}{d}'
    else:
        return f'{n}/{d} {k}' #{sp}{k}{sp}/{sp}{d}'

def ratlin_basic(a: Rational, b: Rational, k: str) -> str:
    if b == 0:
        return rational(a)

    bk = rational(Fraction(abs(b.numerator), abs(b.denominator)), k)
    if a == 0:
        return '-' + bk if b < 0 else bk

    if a.numerator < 0 and b.numerator >= 0:
        return bk + ' - ' + rational(-a)

    return rational(a) + (' + ' if b.numerator >= 0 else ' - ') + bk

def factorout(o: Fraction, a: Rational, b: Rational, k: str) -> str:
    inner = ratlin_basic(a / o, b / o, k)
    if o.numerator == 1 and o.denominator == 1:
        return inner
    if '+' in inner or '-' in inner:
        inner = '(' + inner + ')'
    else:
        inner = ' ' + inner + ' '
    return rational(o, inner)

def ratlinstring(a: Rational, b: Rational, k: str) -> str:
    # Format a + b × k in the most economical way.
    # Try taking out various factors:
    an, ad = abs(a.numerator), a.denominator
    bn, bd = abs(b.numerator), b.denominator
    result = None
    for on in gcd(an, bn), 1, min(an, bn), max(an, bn), lcm(an, bn):
        for od in gcd(ad, bd), 1, min(ad, bd), max(ad, bd), lcm(ad, bd):
            if on != 0 and od != 0:
                for son in -on, on:
                    attempt = factorout(Fraction(son, od), a, b, k)
                    if result == None or len(attempt) < len(result):
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
