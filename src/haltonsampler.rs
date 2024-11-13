use crate::primes::{PRIMES, PRIME_TABLE_SIZE};
use crate::randomizer::DigitPermutation;
use crate::sampler::{RandomStrategy, Sampler, ONE_MINUS_EPSILON};

pub struct HaltonSampler {
    pub index: usize,
    dim: usize,
    strategy: RandomStrategy,
    permuters: Option<Vec<DigitPermutation>>,
}

impl<Real> Sampler<Real> for HaltonSampler
where
    Real: num_traits::Float,
{
    fn get1d(&self) -> Real {
        self.sample_dimension(self.index, self.dim)
    }

    fn get2d(&self) -> [Real; 2] {
        let v1 = self.sample_dimension(self.index, self.dim);
        let v2 = self.sample_dimension(self.index, self.dim + 1);
        [v1, v2]
    }
}

impl HaltonSampler {
    pub fn new(a: usize, dim: usize) -> Self {
        HaltonSampler {
            index: a,
            dim,
            strategy: RandomStrategy::None,
            permuters: None,
        }
    }

    pub fn new_randomized(a: usize, dim: usize, strategy: RandomStrategy) -> Self {
        match strategy {
            RandomStrategy::PermuteDigits => {
                let mut perms = Vec::<DigitPermutation>::new();
                for p in PRIMES.iter().take(PRIME_TABLE_SIZE) {
                    perms.push(DigitPermutation::new(*p as i32, 0));
                }

                HaltonSampler {
                    index: a,
                    dim,
                    strategy,
                    permuters: Some(perms),
                }
            }

            _ => HaltonSampler::new(a, dim),
        }
    }

    fn sample_dimension<Real>(&self, a: usize, dim: usize) -> Real
    where
        Real: num_traits::Float,
    {
        match self.strategy {
            RandomStrategy::PermuteDigits => {
                if let Some(r) = &self.permuters {
                    scramble_radical_inverse(a, dim, &r[dim])
                } else {
                    todo!()
                }
            }
            _ => radical_inverse(a, dim),
        }
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

pub fn scramble_radical_inverse<Real>(
    mut a: usize,
    base_index: usize,
    perm: &DigitPermutation,
) -> Real
where
    Real: num_traits::Float,
{
    assert!(base_index < PRIME_TABLE_SIZE);
    let base = PRIMES[base_index] as usize;
    let inv_base = (Real::one()) / (Real::from(base).unwrap());
    let mut inv_base_m = Real::one();
    //reversed digits:
    let mut rev_digits: usize = 0;
    let mut d_i = 0;
    while Real::from(1 - (base - 1)).unwrap() * inv_base_m < Real::one() {
        let next: usize = a / base;
        // least significant digit
        let digit = (a - next * base) as i32;
        rev_digits = rev_digits * base + perm.permute(d_i, digit) as usize;
        inv_base_m = inv_base_m * inv_base;
        d_i += 1;
        a = next;
    }
    // can be expressed as (d_1*b^(m-1) + d_2*b^(m-2) ... + d_m*b^0 )/b^(m)
    let inv = Real::from(rev_digits).unwrap() * inv_base_m;
    Real::min(inv, Real::from(ONE_MINUS_EPSILON).unwrap())
}

#[test]
fn test_halton_sample() {
    let mut s = HaltonSampler::new_randomized(1, 4, RandomStrategy::PermuteDigits);
    for d in 0..4 {
        println!("---------------------------------------------");
        s.dim = d;
        for _ in 0..16 {
            let v: f32 = s.get1d();
            let p: [f32; 2] = s.get2d();
            s.index += 1;
            print!("{:.2} | ", v);
            println!("{:.2}-{:.2}", p[0], p[1]);
        }
    }
}
