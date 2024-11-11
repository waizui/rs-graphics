// largest float number less than 1
#[cfg(target_pointer_width = "64")]
pub const ONE_MINUS_EPSILON: f64 = 1.0 - f64::EPSILON;
#[cfg(target_pointer_width = "32")]
pub const ONE_MINUS_EPSILON: f32 = 1.0 - f32::EPSILON;

pub trait Sampler<Real>
where
    Real: num_traits::Float,
{
    fn get1d(&self) -> Real;
    fn get2d(&self) -> [Real; 2];
}
