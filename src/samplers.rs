use crate::sobolmatrices::{SOBOL_DIMENSIONS, SOBOL_MATRICES32, SOBOL_MATRIX_SIZE};

pub trait Sampler<Real>
where
    Real: num_traits::Float,
{
    fn get1d(&self) -> Real;
    fn get2d(&self) -> [Real; 2];
}

fn reverse_bit_32(mut n: u32) -> u32 {
    n = (n << 16) | (n >> 16);
    n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
    n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
    n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
    n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
    n
}

/// elements c are integers represent column of generator matrix
fn mul_generator(c: &[u32], d: usize) -> u32 {
    let mut v = 0;
    let mut i = 0;
    let mut a = d;
    while a != 0 {
        if a & 1 == 1 {
            v ^= c[i];
        }

        a >>= 1;
        i += 1;
    }
    v
}

pub struct SobolSampler {
    a: usize,
    dim: usize,
}

fn sobol_sample<Real>(a: usize, dim: usize) -> Real
where
    Real: num_traits::Float,
{
    assert!(dim < SOBOL_DIMENSIONS);
    assert!(a < (1 << SOBOL_MATRIX_SIZE));
    // not supporting 64bit yet
    let v = mul_generator(&SOBOL_MATRICES32, dim);
    todo!()
}

impl<Real> Sampler<Real> for SobolSampler
where
    Real: num_traits::Float,
{
    fn get1d(&self) -> Real {
        sobol_sample::<Real>(self.a, self.dim)
    }

    fn get2d(&self) -> [Real; 2] {
        let v1 = sobol_sample::<Real>(self.a, self.dim);
        let v2 = sobol_sample::<Real>(self.a, self.dim + 1);
        [v1, v2]
    }
}
