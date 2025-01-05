
const LIMIT: usize = 0xffffffff;
// n is stored at index (n-3)/2.  Index i represents 2i + 3
const BITS: usize = (LIMIT - 3) / 2;

//static mut SSIEVE: [bool; BITS]; // = [false; BITS];
use std::io::Write;

pub fn main() {
    let mut sieve = Box::new([false; BITS]);// unsafe { &mut SSIEVE };
    println!("2");
    let mut binary = std::io::BufWriter::new(std::fs::File::create("primes.bin").unwrap());
    binary.write(&2u32.to_le_bytes()).unwrap();
    let mut n = 0;
    loop {
        let p = 2*n + 3;
        if p * p > LIMIT {
            break;
        }
        if sieve[n] {
            n = n + 1;
            continue;
        }
        n = n + 1;
        println!("{}", p);
        binary.write(&(p as u32).to_le_bytes()).unwrap();
        for i in ((p*p - 3)/2 .. BITS).step_by(p) {
            sieve[i] = true;
        }
    }
    for i in n .. BITS {
        if i & 65535 == 0 {
            println!("{}", i);
        }
        if !sieve[i] {
            let p = 2*i + 3;
            // println!("{}", p);
            binary.write(&(p as u32).to_le_bytes()).unwrap();
        }
    }
}
