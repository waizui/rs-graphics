use crate::primes::{PRIMES, PRIME_TABLE_SIZE};
use crate::sampler::{Sampler, ONE_MINUS_EPSILON};

pub struct HaltonSampler {
    pub a: usize,
    pub dim: usize,
}

impl<Real> Sampler<Real> for HaltonSampler
where
    Real: num_traits::Float,
{
    fn get1d(&self) -> Real {
        radical_inverse(self.a, self.dim)
    }

    fn get2d(&self) -> [Real; 2] {
        let v1 = radical_inverse(self.a, self.dim);
        let v2 = radical_inverse(self.a, self.dim + 1);
        [v1, v2]
    }
}

pub fn radical_inverse<Real>(mut a: usize, base_index: usize) -> Real
where
    Real: num_traits::Float,
{
    assert!(base_index < PRIME_TABLE_SIZE);
    let base = PRIMES[base_index] as usize;
    let inv_base = (Real::one()) / (Real::from(base).unwrap());
    let mut inv_base_m = Real::one();
    //reversed digits:
    let mut rev_digits: usize = 0;
    while a != 0 {
        let next: usize = a / base;
        // least significant digit
        let digit: usize = a - next * base;
        rev_digits = rev_digits * base + digit;
        inv_base_m = inv_base_m * inv_base;
        a = next;
    }
    // can be expressed as (d_1*b^(m-1) + d_2*b^(m-2) ... + d_m*b^0 )/b^(m)
    let inv = Real::from(rev_digits).unwrap() * inv_base_m;
    Real::min(inv, Real::from(ONE_MINUS_EPSILON).unwrap())
}

pub fn reverse_bit_32(mut n: u32) -> u32 {
    n = (n << 16) | (n >> 16);
    n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
    n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
    n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
    n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
    n
}

#[test]
fn test_halton_sample() {
    let mut s = HaltonSampler { a: 1, dim: 1 };
    for d in 0..4 {
        println!("---------------------------------------------");
        s.dim = d;
        for _ in 0..16 {
            let v: f32 = s.get1d();
            let p: [f32; 2] = s.get2d();
            s.a += 1;
            print!("{:.2} | ", v);
            println!("{:.2}-{:.2}", p[0], p[1]);
        }
    }
}
