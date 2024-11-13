use std::vec::Vec;


pub struct DigitPermutation {
    ndigits: i32,
    base: i32,
    permutations: Vec<u16>,
}

impl DigitPermutation {
    pub fn new(base: i32, seed: i32) -> Self {
        // permutations storing u16
        assert!(base < 65536);
        let mut ndigits: i32 = 0;

        // number of digits needed for base
        let inv_base = 1. / (base as f32);
        let mut inv_base_m: f32 = 1.;
        loop {
            if 1. - ((base as f32) - 1.) * inv_base_m < 1. {
                // is the least significant digit
                break;
            }
            ndigits += 1;
            inv_base_m *= inv_base;
        }
        let size: usize = (ndigits * base) as usize;
        let mut permutations: Vec<u16> = vec![0u16; size];
        // compute random permutations for all digits
        let mut d_i: i32 = 0;
        while d_i < ndigits {
            d_i += 1;
            let dseed = hash(base, d_i, seed);
            let mut d_v = 0;
            while d_v < base {
                let i = (d_i * base + d_v) as usize;
                permutations[i] = permutation_elements(d_v as u32, base as u32, dseed) as u16;
                d_v += 1;
            }
        }

        DigitPermutation {
            ndigits,
            base,
            permutations,
        }
    }

    pub fn permute(&self, digit_i: i32, digit_v: i32) -> i32 {
        assert!(digit_i < self.ndigits);
        assert!(digit_v < self.base);
        let i = (digit_i * self.base + digit_v) as usize;
        self.permutations[i].into()
    }
}

/// hash function from: https://graphics.pixar.com/library/MultiJitteredSampling/paper.pdf, Andrew Kensler
/// i is value, l is base , p is seed
pub fn permutation_elements(mut i: u32, l: u32, p: u32) -> i32 {
    // we want a num in range [0,l-1]
    let mut w = l - 1;
    // filling 1 to the right higest 1, efficient:  00100100 -> 00111111
    w |= w >> 1;
    w |= w >> 2;
    w |= w >> 4;
    w |= w >> 8;
    w |= w >> 16;

    loop {
        i ^= p;
        i *= 0xe170893d;
        i ^= p >> 16;
        i ^= (i & w) >> 4;
        i ^= p >> 8;
        i *= 0x0929eb3f;
        i ^= p >> 23;
        i ^= (i & w) >> 1;
        i *= 1 | p >> 27;
        i *= 0x6935fa69;
        i ^= (i & w) >> 11;
        i *= 0x74dcb303;
        i ^= (i & w) >> 2;
        i *= 0x9e501cc3;
        i ^= (i & w) >> 2;
        i *= 0xc860a3df;
        i &= w;
        i ^= i >> 5;

        if i >= l {
            break;
        }
    }

    ((i + p) % l) as i32
}

// hash three values
fn hash(v1: i32, v2: i32, v3: i32) -> u32 {
    let mut buf = [0u8; 12];
    let mut i: usize = 0;
    for v in [v1, v2, v3] {
        buf[i] = (v >> 24) as u8;
        buf[i + 1] = (v >> 16) as u8;
        buf[i + 2] = (v >> 8) as u8;
        buf[i + 3] = v as u8;
        i += 4;
    }

    murmur_hash_64a(&buf, 0) as u32
}

// https://github.com/explosion/murmurhash/blob/master/murmurhash/MurmurHash2.cpp
pub fn murmur_hash_64a(key: &[u8], seed: u64) -> u64 {
    const M: u64 = 0xc6a4a7935bd1e995;
    const R: i32 = 47;

    let len = key.len();
    let mut h = seed ^ (len as u64).wrapping_mul(M);

    // chunks of 8 bytes
    for chunk in 0..len / 8 {
        let offset = chunk * 8;
        // compose k from 8 bytes
        let mut k = u64::from(key[offset])
            | (u64::from(key[offset + 1]) << 8)
            | (u64::from(key[offset + 2]) << 16)
            | (u64::from(key[offset + 3]) << 24)
            | (u64::from(key[offset + 4]) << 32)
            | (u64::from(key[offset + 5]) << 40)
            | (u64::from(key[offset + 6]) << 48)
            | (u64::from(key[offset + 7]) << 56);

        k = k.wrapping_mul(M);
        k ^= k >> R;
        k = k.wrapping_mul(M);

        h ^= k;
        h = h.wrapping_mul(M);
    }

    // remaining bits
    let r = len & 7;
    if r > 0 {
        let offset = len - r;
        match r {
            7 => h ^= u64::from(key[offset + 6]) << 48,
            6 => h ^= u64::from(key[offset + 5]) << 40,
            5 => h ^= u64::from(key[offset + 4]) << 32,
            4 => h ^= u64::from(key[offset + 3]) << 24,
            3 => h ^= u64::from(key[offset + 2]) << 16,
            2 => h ^= u64::from(key[offset + 1]) << 8,
            1 => h ^= u64::from(key[offset]),
            _ => unreachable!(),
        }
        h = h.wrapping_mul(M);
    }

    h ^= h >> R;
    h = h.wrapping_mul(M);
    h ^= h >> R;

    h
}

#[test]
fn test_murmur_hash() {
    let test_data = b"this is a random text";
    let seed = 0x12345678;
    let hash = murmur_hash_64a(test_data, seed);
    assert_ne!(hash, 0); // sanity check
}
