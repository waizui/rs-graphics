type Real = f32;
const PI: Real = std::f32::consts::PI;

pub trait Image {
    fn new(w: usize, h: usize) -> Self;
    fn w(&self) -> usize;
    fn h(&self) -> usize;
    fn get(&self, uv: &[Real; 2]) -> [Real; 3];
    fn set(&mut self, uv: &[Real; 2], color: &[Real; 3]);

    fn get_p(&self, uv: &[usize; 2]) -> [Real; 3] {
        let p = pixel2tex(uv[0], uv[1], self.w(), self.h());
        self.get(&p)
    }
    fn set_p(&mut self, uv: &[usize; 2], color: &[Real; 3]) {
        let p = pixel2tex(uv[0], uv[1], self.w(), self.h());
        self.set(&p, color)
    }
}

pub fn roundi(f: Real) -> i32 {
    (f + 0.5) as i32
}

pub fn tex2pixel(t: Real, ext: usize) -> usize {
    roundi(t * ((ext - 1) as Real) - 0.5) as usize
}

pub fn pixel2tex(x: usize, y: usize, w: usize, h: usize) -> [Real; 2] {
    let u = x as Real / w as Real;
    let v = y as Real / h as Real;
    [u, v]
}

pub fn unitsphere2envmap(d: &[f32; 3]) -> [Real; 2] {
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

// https://github.com/mmp/pbrt-v4/blob/1ae72cfa7344e79a7815a21ed3da746cdccee59b/src/pbrt/util/math.cpp#L292
pub fn envmap2unitsphere(p: &[Real; 2]) -> [Real; 3] {
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

fn calc_grayscale<T>(img: &T) -> T
where
    T: Image,
{
    let w = img.w();
    let h = img.h();
    let mut grayscale = T::new(w, h);

    for i_h in 0..h {
        for i_w in 0..w {
            let p = [i_w, i_h];
            let c = img.get_p(&p);
            let gray = (0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2]).clamp(0., 1.);
            grayscale.set_p(&p, &[gray, gray, gray]);
        }
    }

    grayscale
}

fn calc_integral_over_grayscale<T>(gray_img: &T) -> Real
where
    T: Image,
{
    let w = gray_img.w();
    let h = gray_img.h();

    let mut sum = 0.;

    for i_h in 0..h {
        for i_w in 0..w {
            let c = gray_img.get_p(&[i_w, i_h])[0];
            sum += c;
        }
    }

    sum /= (w * h) as Real;
    sum
}

pub fn calc_avg_col_row<T>(grayscale: &T) -> [T; 2]
where
    T: Image,
{
    let w = grayscale.w();
    let h = grayscale.h();
    let mut col_avg_map = T::new(w, 1);

    for p_w in 0..w {
        let mut c_avg = 0.;
        for p_h in 0..h {
            let gray = grayscale.get_p(&[p_w, p_h])[0];
            c_avg += gray;
        }
        c_avg /= h as Real;
        col_avg_map.set_p(&[p_w, 1], &[c_avg; 3]);
    }

    let mut row_avg_map = T::new(1, h);
    for p_h in 0..h {
        let mut r_avg = 0.;
        for p_w in 0..w {
            let gray = grayscale.get_p(&[p_w, p_h])[0];
            r_avg += gray;
        }
        r_avg /= w as Real;
        row_avg_map.set_p(&[1, p_h], &[r_avg; 3]);
    }

    [col_avg_map, row_avg_map]
}

/// use tow channels to store cdf^-1 on w and h
pub fn calc_inverse_cdf_map<T>(envmap: &T) -> T
where
    T: Image,
{
    let grayscale = calc_grayscale(envmap);
    let avgs = calc_avg_col_row(&grayscale);
    let col_avg = &avgs[0];
    let row_avg = &avgs[1];

    let w = col_avg.w();
    let h = row_avg.h();

    let mut map = T::new(w, h);

    for p_w in 0..w {
        for p_h in 0..h {
            let mut sum = 0.;

            let mut i_x = 0;
            let u = p_w as Real / w as Real;
            for i_w in 0..w {
                // p(w|h) = p(w,h)/p(h) = gray(w,h) / col_avg(w)
                let c_avg = col_avg.get_p(&[i_w, 0])[0]; // col avg of w-th column
                sum += grayscale.get_p(&[i_w, p_h])[0] / (c_avg * w as Real);
                if sum >= u {
                    i_x = i_w;
                    break;
                }
            }

            sum = 0.;
            let mut i_y = 0;
            let v = p_h as Real / w as Real;
            for i_h in 0..h {
                let r_avg = row_avg.get_p(&[0, p_h])[0]; // row avg of h-th row
                sum += grayscale.get_p(&[p_w, i_h])[0] / (r_avg * h as Real);
                if sum >= v {
                    i_y = i_h;
                    break;
                }
            }

            map.set_p(
                &[p_w, p_h],
                &[i_x as Real / w as Real, i_y as Real / h as Real, 0.],
            );
        }
    }

    map
}

#[test]
fn test_sphere_mapping() {
    let mut dir = [0.; 3];

    for i in 0..100 {
        for j in 0..3 {
            dir[j] = (i + 1) as Real;
        }

        del_geo_core::vec3::normalize(&mut dir);

        let uv = unitsphere2envmap(&dir);
        let udir = envmap2unitsphere(&uv);

        assert!((dir[0] - udir[0]).abs() < 0.001);
        assert!((dir[1] - udir[1]).abs() < 0.001);
        assert!((dir[2] - udir[2]).abs() < 0.001);
    }
}
