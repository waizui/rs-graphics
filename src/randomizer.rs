use std::vec::Vec;

pub struct DigitPermutation {
    ndigits: u32,
    base: u32,
    permutations: Vec<u16>,
}

impl DigitPermutation {
    pub fn new(base: u32, seed: u32) -> Self {
        // permutations storing u16
        assert!(base < 65536);
        let mut ndigis: u32 = 0;

        // number of digits needed for base
        let inv_base = 1. / (base as f32);
        let mut inv_base_m: f32 = 1.;
        loop {
            if 1. - ((base as f32) - 1.) * inv_base_m < 1. {
                // is the least significant digit
                break;
            }
            ndigis += 1;
            inv_base_m *= inv_base;
        }
        let size: usize = (ndigis * base) as usize;
        let mut permutations: Vec<u16> = vec![0u16; size];
        // compute random permutations for all digits
        let mut d_i: u32 = 0;
        while d_i < ndigis {
            d_i += 1;
            let dseed = hash(base, d_i, seed);
            let mut d_v = 0;
            while d_v < base {
                let i = (d_i * base + d_v) as usize;
                permutations[i] = permutation_elements(d_v, base, dseed) as u16;
                d_v += 1;
            }
        }

        DigitPermutation {
            ndigits: ndigis,
            base: base,
            permutations,
        }
    }

    pub fn permute(&self, digit_i: u32, digit_v: u32) -> i32 {
        assert!(digit_i < self.ndigits);
        assert!(digit_v < self.base);
        let i = (digit_i * self.base + digit_v) as usize;
        self.permutations[i].into()
    }
}

/// hash function from: https://graphics.pixar.com/library/MultiJitteredSampling/paper.pdf, Andrew Kensler
/// i is value, l is base , p is seed
pub fn permutation_elements(mut i: u32, mut l: u32, p: u32) -> i32 {
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

fn hash(mut v: u32, mut i: u32, p: u32) -> u32 {
    //TODO:impl
    todo!()
}
