use crate::sampler::Sampler;
use crate::sobolmatrices::{SOBOL_DIMENSIONS, SOBOL_MATRICES32, SOBOL_MATRIX_SIZE};

pub struct SobolSampler {
    pub a: usize,
    pub dim: usize,
}

impl<Real> Sampler<Real> for SobolSampler
where
    Real: num_traits::Float,
{
    fn get1d(&self) -> Real {
        sobol_sample(self.a, self.dim)
    }

    fn get2d(&self) -> [Real; 2] {
        let v1 = sobol_sample(self.a, self.dim);
        let v2 = sobol_sample(self.a, self.dim + 1);
        [v1, v2]
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
