type Real = f32;
const PI: Real = std::f32::consts::PI;

pub trait Image {
    fn w(&self) -> Real;
    fn h(&self) -> Real;
    fn get(&self, uv: &[Real; 2]) -> [Real; 3];
    fn set(&self, uv: &[Real; 2], color: &[Real; 3]);
}

pub fn gen_inverse_cdf_map<T>(img: &T) -> T
where
    T: Image,
{
    let w = img.w();
    let h = img.h();
    todo!()
}

