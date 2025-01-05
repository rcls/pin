
const FACTORS : [u64; 7] = [3, 5, 17, 257, 65537, 641, 6700417];

fn double(a: u64, p: u64) -> u64 {
    if a & (1 << 63) != 0 {
        (a ^ 1 << 63) << 1 ^ p
    }
    else {
        a << 1
    }
}

fn mult(mut a: u64, mut b: u64, p: u64) -> u64 {
    let mut res = 0;
    while a != 0 {
        if a & 1 != 0 {
            res ^= b;
        }
        a >>= 1;
        b = double(b, p);
    }
    res
}

fn square(a: u64, p: u64) -> u64 {
    mult(a, a, p)
}

fn power(mut a: u64, mut e: u64, p: u64) -> u64 {
    let mut res = 1;
    while e != 0 {
        if e & 1 != 0 {
            res = mult(res, a, p);
        }
        e >>= 1;
        a = square(a, p);
    }
    res
}

fn is_cyclic(p: u64) -> bool {
    const O: u64 = 0xffffffffffffffff;
    if power(2, O, p) != 1 {
        return false;
    }
    for f in FACTORS {
        assert!(O % f == 0);
        if power(2, O / f, p) == 1 {
            return false;
        }
    }
    true
}

fn scan(bits: usize, base: u64, end: u64) {
    if bits == 0 {
        if is_cyclic(base) {
            println!("{:#018x}", base);
        }
        return;
    }
    for i in 1 .. end {
        scan(bits - 1, base | 1 << i, i)
    }
}

pub fn main() {
    scan(1, 1, 64);
    scan(3, 1, 64);
}
