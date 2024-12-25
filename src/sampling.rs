type Real = f32;
const PI: Real = std::f32::consts::PI;


/// returen  [-1,1]
pub fn sample_uni_disk_concentric(p: &[Real; 2]) -> [Real; 2] {
    let x = p[0] * 2. - 1.;
    let y = p[1] * 2. - 1.;
    if x == 0. && y == 0. {
        return [0.; 2];
    }

    let (theta, r) = {
        if x.abs() > y.abs() {
            (x, PI / 4. * (y / x))
        } else {
            (x, PI / 4. * (x / y))
        }
    };

    [r * theta.cos(), r * theta.sin()]
}
