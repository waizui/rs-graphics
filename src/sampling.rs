use del_geo_core::vec3::Vec3;

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

/// d: point on unit sphere surface
/// return u,v in [0,1]^2
pub fn unitsphere2square(d: &[f32; 3]) -> [Real; 2] {
    let x = d[0].abs();
    let y = d[1].abs();
    let z = d[2].abs();
    let r = (1. - z).sqrt();
    let phi = y.atan2(x);
    let phi = phi * std::f32::consts::FRAC_2_PI;
    let v = phi * r;
    let u = r - v;
    let (u, v) = if d[2] < 0. { (1. - v, 1. - u) } else { (u, v) };
    let u = u.copysign(d[0]);
    let v = v.copysign(d[1]);
    [u * 0.5 + 0.5, v * 0.5 + 0.5]
}

/// p: point on unit square
/// return [x,y,z] on unit sphere surface
// https://github.com/mmp/pbrt-v4/blob/1ae72cfa7344e79a7815a21ed3da746cdccee59b/src/pbrt/util/math.cpp#L292
pub fn square2unitsphere(p: &[Real; 2]) -> [Real; 3] {
    // map [0,1] to [-1,1]
    let u = 2. * p[0] - 1.;
    let v = 2. * p[1] - 1.;
    let up = u.abs();
    let vp = v.abs();
    let sd = 1. - (up + vp);
    let d = sd.abs();
    let r = 1. - d;
    let phi = (if r == 0. { 1. } else { (vp - up) / r + 1. }) * PI / 4.;
    let z = (1. - r * r).copysign(sd);
    let cosphi = phi.cos().copysign(u);
    let sinphi = phi.sin().copysign(v);

    let x = cosphi * r * (2. - r * r).sqrt();
    let y = sinphi * r * (2. - r * r).sqrt();

    [x, y, z]
}
