#![feature(cfg_boolean_literals)]

// use std::collections::hash_set::HashSet;

pub fn print_table(f : impl Fn(usize, usize) -> usize,
                   cols: impl Clone + Iterator<Item=usize>,
                   rows: impl Iterator<Item=usize>) {
    print!("\x1b[4m *\u{2502}");
    let mut j = cols.clone();
    print!("{:2}", j.next().unwrap());
    for b in j {
        print!("{b:3}");
    }
    println!("\x1b[m");
    for a in rows {
        let mut j = cols.clone();
        print!("{a:2}\u{2502}{:2}", f(a, j.next().unwrap()));
        for b in j {
            print!("{:3}", f(a, b));
        }
        println!();
    }
}

pub trait Laver {
    fn l(&self, a: usize, b: usize) -> usize;
    fn o(&self, a: usize, b: usize) -> usize {
        let n = self.len();
        let b1 = if b == n {1} else {b+1};
        let o1 = self.l(a, b1);
        if o1 == 1 {n} else {o1 - 1}
    }
    fn len(&self) -> usize;

    fn print_table(&self) {
        let n = self.len();
        let f = |x| if x == 0 {n} else {x};
        print_table(|a,b| self.l(f(a),f(b)) % n,
                    (1..=n).map(|x| x % n),
                    (1..=n).map(|y| y % n));
    }

    fn subtable(&self, scale: usize) -> Subtable<&Self> {
        Subtable::new(scale, self)
    }
}

impl<A: Laver + ?Sized> Laver for &A {
    fn len(&self) -> usize { (*self).len() }
    fn l(&self, a: usize, b: usize) -> usize { (*self).l(a, b) }
}

pub fn check_eq<T: Laver, U: Laver>(t: &T, u: &U) {
    let n = t.len();
    assert_eq!(n, u.len());
    for a in 1 ..= n {
        for b in 1 ..= n {
            assert_eq!(t.l(a,b), u.l(a,b), "{a} {b}");
        }
    }
}

pub struct Table {
    n : usize,
    cycle: Vec<Vec<usize>>,
}

impl Table {
    pub fn new(n: usize) -> Table {
        let mut cycle = Vec::new();
        cycle.resize_with(n + 1, || Vec::new());
        let mut r = Table{n, cycle};
        let last = &mut r.cycle[n];
        last.push(n);
        for i in 1..n {
            last.push(i);
        }

        for a in (1..n).rev() {
            r.cycle[a].push(n);
            let mut next = a+1;
            while next != n {
                r.cycle[a].push(next);
                next = r.l(next, a+1);
            }
        }
        r
    }
    pub fn print_cycles(&self) {
        for a in 1..self.n {
            println!("{:2}: {:x?}", a, &self.cycle[a]);
        }
    }
}

impl Laver for Table {
    fn len(&self) -> usize { self.n }
    fn l(&self, a: usize, b: usize) -> usize {
        assert!(a > 0);
        assert!(b > 0);
        assert!(a <= self.n);
        assert!(b <= self.n);
        let aa = &self.cycle[a];
        aa[b % aa.len()]
    }
}

pub struct Periods {
    n : usize,
    periods: Vec<(usize, usize)>,
}

impl Periods {
    pub fn new(n: usize) -> Periods {
        assert_eq!(n & (n-1), 0);
        let mut periods = Vec::new();
        periods.resize(n+1, (0, 0));
        let mut r = Periods{n, periods};
        if n == 1 {
            return r;
        }
        r.periods[n] = (n, n/2);
        r.periods[n-1] = (1, 0);
        'create: for a in (1..n-1).rev() {
            // Change leading 0 to 1...
            let l1 = (r.n - 1 - a).leading_zeros();
            let bit = 1 << (usize::BITS - 1 - l1);
            let aa = a | bit;
            assert!(aa > a, "{a} {aa}");
            assert!(aa < n, "{a} {aa}");
            let (pp, _) = r.periods[aa];

            let mut v = a + 1;
            for i in 0 .. pp {
                if v > aa {
                    r.periods[a] = (pp, i);
                    continue 'create;
                }
                v = r.l(v, a + 1);
            }
            assert!(v == aa + 1);
            r.periods[a] = (2 * pp, pp);
        }
        r
    }

    pub fn period(&self, a: usize) -> usize {
        self.periods[a].0
    }
}

impl Laver for Periods {
    fn len(&self) -> usize { self.n }
    // Specific to power of 2.
    fn l(&self, mut a: usize, mut b: usize) -> usize {
        if a == self.n {
            return b;
        }
        // println!("< {} {}", a, b);
        let mut sub = 0;
        loop {
            if a == self.n - 1 {
                break;
            }
            let (p, theta) = self.periods[a];
            b = b & (p-1);
            if b == 0 {
                break;
            }
            let l1 = (self.n - 1 - a).leading_zeros();
            let bit = 1 << (usize::BITS - 1 - l1);
            // println!(". {} {} {} {}", a, b, bit, self.jumps[a]);
            if b <= theta {
                sub += bit;
            }
            a += bit;
        }
        // println!("> {sub}");
        self.n - sub
    }
}

pub struct Subtable<T> {
    n: usize,
    scale: usize,
    source: T
}

impl<T: Laver> Subtable<T> {
    pub fn new(scale: usize, source: T)  -> Self {
        let n = source.len() / scale;
        assert_eq!(n * scale, source.len());
        Subtable{n, scale, source}
    }
}

impl<T: Laver> Laver for Subtable<T> {
    fn len(&self) -> usize { self.n }
    fn l(&self, a: usize, b: usize) -> usize {
        self.source.l(a * self.scale, b * self.scale) / self.scale
    }
}

pub fn main() {
    println!("{:?}", std::env::args());

    for arg in std::env::args().skip(1) {
        let n = usize::from_str_radix(&arg, 10).unwrap();

        //let laver = Table::new(n);
        //laver.print_cycles();
        //if n <= 1024 {
        //    laver.print_table();
        //}
        let l = Periods::new(n);
        //check_eq(&laver, &l);
        for (i, pj) in l.periods[1..].iter().enumerate() {
            println!("{i} {pj:?}");
        }
        if n <= 1024 {
            l.print_table();
        }
        //for i in 3..30 {
        //    println!("{i} {}", Table::new(1 << i).cycle[5].len());
        //}
        println!();
        println!("Matches...");
    }

    //check_eq(&Periods::new(4), &Periods::new(32).subtable(8));
    //check_eq(&Periods::new(64), &Periods::new(1024).subtable(16));
    //check_eq(&Periods::new(128), &Periods::new(2048).subtable(16));
}

#[cfg(false)]
fn check_final() { // Passes
    let n = 1 << 25;
    let l = Table::new(n);
    for i in 1 .. n/2 {
        let bit = !i & (i+1);
        assert!(bit & (bit - 1) == 0);
        assert!((i+1) & bit != 0);
        assert!((i+1) & (bit-1) == 0);
        assert!(l.l(i, n-1) >= n - bit, "{n} {i}");
    }
}

#[cfg(false)]
fn check_jump() {
    for i in 2..10 {
        println!("{i}");
        let n = 1 << i;
        let hn = 1 << (i-1);
        let l = Periods::new(n);
        for a in 1 .. hn {
            println!("{i} {a}");
            if a & a - 1 == 0 {
                continue;
            }
            let mut b = l.periods[a];
            assert!(b & b - 1 == 0);
            if b <= 2 {
                continue;
            }
            while l.l(a, b-1) >= hn {
                b = b / 2;
                assert!(b & b - 1 == 0);
                assert!(b != 0);
            }
            println!("{i} {a} {b}");
            assert!(l.l(a, b) >= hn, "{} {}", l.l(a,b-1), l.l(a,b));
        }
    }
}
