use crate::randomizer::{hash, permutation_elements};
use crate::sampler::{RandomStrategy, Sampler};
use crate::sobolmatrices::{SOBOL_DIMENSIONS, SOBOL_MATRICES32, SOBOL_MATRIX_SIZE};

pub struct SobolSampler {
    pub index: usize,
    pub dim: usize,
    strategy: RandomStrategy,
}

impl<Real> Sampler<Real> for SobolSampler
where
    Real: num_traits::Float,
{
    fn set_i(&mut self, i: usize) {
        self.index = i;
    }

    fn set_dim(&mut self, dim: usize) {
        self.dim = dim;
    }

    fn get1d(&mut self) -> Real {
        let dim = self.dim;
        self.dim += 1;
        // parameters to hash could be adjust for meanful
        let hash = hash(self.index as i32, dim as i32, 0);
        // 1024 is the max sample count
        let index = permutation_elements(self.index as u32, 1024, hash as u32) as usize;
        sample_dimension(&self.strategy, index, dim, hash as usize)
    }

    fn get2d(&mut self) -> [Real; 2] {
        let dim = self.dim;
        self.dim += 2;

        let hash = hash(self.index as i32, dim as i32, 0);
        let index = permutation_elements(self.index as u32, 1024, hash as u32) as usize;
        let v1 = sample_dimension(&self.strategy, index, dim, hash as usize);
        let v2 = sample_dimension(&self.strategy, index, dim + 1, (hash >> 32) as usize);
        [v1, v2]
    }

    fn restore(&mut self) {
        self.dim = 0;
        self.index = 1;
    }
}

impl Default for SobolSampler {
    fn default() -> Self {
        SobolSampler {
            index: 1,
            dim: 0,
            strategy: RandomStrategy::None,
        }
    }
}

impl SobolSampler {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn new_randomized(strategy: RandomStrategy) -> Self {
        SobolSampler {
            index: 1,
            dim: 0,
            strategy,
        }
    }
}

fn sample_dimension<Real>(strategy: &RandomStrategy, a: usize, dim: usize, hash: usize) -> Real
where
    Real: num_traits::Float,
{
    match strategy {
        // TODO: explanary comment
        RandomStrategy::PermuteDigits => sobol_sample(a, dim, |x| hash ^ x),
        _ => sobol_sample(a, dim, |x| x),
    }
}

/// r: randomizer
pub fn sobol_sample<Real, F>(mut a: usize, dim: usize, r: F) -> Real
where
    Real: num_traits::Float,
    F: Fn(usize) -> usize,
{
    assert!(dim < SOBOL_DIMENSIONS);
    assert!(std::mem::size_of::<usize>() == 4 || a < (1 << SOBOL_MATRIX_SIZE));

    let mut v = 0;
    let mut i = dim * SOBOL_MATRIX_SIZE;
    // can be expressed as: v = d_1(a)*c_1 + d_2(a)*c_2 ...d_32(a)*c_32, where d_i(a) is the i-th digit of a,
    // c_i represents column of generator matrix.
    // eg. 4 = 100, matrix = [4,2,1] => [4,2,1][0,0,1]^t = 001
    while a != 0 {
        if a & 1 == 1 {
            // bitwise add
            v ^= SOBOL_MATRICES32[i];
        }

        a >>= 1;
        i += 1;
    }

    // 0x2f800000 = 1^-32f
    let f = (r(v as usize) as f32) * f32::from_bits(0x2f800000);
    Real::from(f).expect("can not convert sobol sample value")
}
