trait Sampler<Real> {
    fn get1D() -> Real;
    fn get2D() -> [Real; 2];
}

trait SobolSampler<Real> {}

impl<T: SobolSampler<Real>, Real> Sampler<Real> for T {
    fn get1D() -> Real {
        todo!()
    }

    fn get2D() -> [Real; 2] {
        todo!()
    }
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
fn mul_generator(c: &[u32], d: u32) -> u32 {
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
