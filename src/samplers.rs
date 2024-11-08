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

pub struct HaltonSampler {}

impl<Real> Sampler<Real> for HaltonSampler
where
    Real: num_traits::Float,
{
    fn get1d(&self) -> Real {
        todo!()
    }

    fn get2d(&self) -> [Real; 2] {
        todo!()
    }
}

pub struct SobolSampler {
    a: usize,
    dim: usize,
}

fn sobol_sample<Real>(mut a: usize, dim: usize) -> Real
where
    Real: num_traits::Float,
{
    assert!(dim < SOBOL_DIMENSIONS);
    assert!(std::mem::size_of::<usize>() == 4 || a < (1 << SOBOL_MATRIX_SIZE));

    let mut v = 0;
    let mut i = dim * SOBOL_MATRIX_SIZE;
    while a != 0 {
        if a & 1 == 1 {
            v ^= SOBOL_MATRICES32[i];
        }

        a >>= 1;
        i += 1;
    }

    let f = (v as f32) * f32::from_bits(0x2f800000);
    Real::from(f).expect("can not convert sobol sample value")
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

#[test]
fn test_sobol_sample() {
    let mut s = SobolSampler { a: 1, dim: 1 };

    for d in 0..4 {
        println!("---------------------------------------------");
        s.dim = d;
        for i in 0..16 {
            let v: f32 = s.get1d();
            let p: [f32; 2] = s.get2d();
            s.a += 1;
            print!("{:.2} | ", v);
            println!("{:.2}-{:.2}", p[0], p[1]);
        }
    }
}
