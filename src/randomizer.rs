use std::vec::Vec;

pub struct Randomizer {
    p: u32,
    l: u32,
    permutations: Vec<u16>,
}

impl Randomizer {
    pub fn new(base: u32, seed: u32) -> Self {
        // permutations storing u16
        assert!(base < 65536);
        let mut ndigis = 0;

        // TODO: calculate permutations

        Randomizer {
            p: seed,
            l: base,
            permutations: Vec::new(),
        }
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
