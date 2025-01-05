
use memmap;

const LIMIT: usize = 0x100000001;

const MAX_BIT: usize = (LIMIT - 3) / 2;
const MAX_BYTE: usize = MAX_BIT / 8;

#[inline]
fn clear_bit(bits: &mut [u8], n: usize) {
    let index = (n-3) / 2;
    bits[index >> 3] &= !(1 << (index & 7));
}

#[inline]
fn get_bit(bits: &[u8], n: usize) -> bool {
    let index = (n-3) / 2;
    bits[index >> 3] & 1 << (index & 7) != 0
}

pub fn main() {
    assert!(MAX_BIT & 7 == 7);
    let file = std::fs::OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(true)
        .open("prime-mask.bin").unwrap();
    file.set_len(MAX_BYTE as u64 + 1).unwrap();
    let mut bits = unsafe { memmap::MmapMut::map_mut(&file).unwrap() };
    let bits = bits.as_mut();
    bits.fill(0xff);

    for n in (3 .. (LIMIT as f64).sqrt() as usize).step_by(2) {
        if get_bit(bits, n) {
            for i in (n * n ..= LIMIT).step_by(2 * n) {
                clear_bit(bits, i);
            }
        }
    }
}
