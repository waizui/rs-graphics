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
        sobol_sample(self.index, dim)
    }

    fn get2d(&mut self) -> [Real; 2] {
        let dim = self.dim;
        self.dim += 2;
        let v1 = sobol_sample(self.index, dim);
        let v2 = sobol_sample(self.index, dim + 1);
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

pub fn sobol_sample<Real>(mut a: usize, dim: usize) -> Real
where
    Real: num_traits::Float,
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
    let f = (v as f32) * f32::from_bits(0x2f800000);
    Real::from(f).expect("can not convert sobol sample value")
}
