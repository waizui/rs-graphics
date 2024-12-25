use del_geo_core::vec3;
use nalgebra::Matrix4;

use crate::{mat4::transform_vector, sampling};

type Real = f32;

#[derive(Debug)]
pub struct Camera {
    pub w: usize,
    pub h: usize,
    pub fov: Real,
    pub n: Real,
    pub f: Real,
    pub pos: [Real; 3],
    pub up: [Real; 3],
    pub forward: [Real; 3],
    pub len_r: Option<Real>,
}

impl Default for Camera {
    fn default() -> Self {
        Camera::new()
    }
}

impl Camera {
    pub fn new() -> Camera {
        let pos = [0.; 3];
        let forward = [0., 0., -1.];
        let n = 0.01;
        let f = 1000.;
        let fov = 45.;
        let up = [0., 1., 0.];
        Camera {
            w: 512,
            h: 512,
            fov,
            n,
            f,
            pos,
            up,
            forward,
            len_r: None,
        }
    }
}

pub fn gen_ray(
    (ix, iy): (usize, usize),
    (dx, dy): (f32, f32),
    img_shape: (usize, usize),
    fov: f32,
    transform_camlcl2world: &[f32; 16],
) -> ([f32; 3], [f32; 3]) {
    assert!(ix < img_shape.0 && iy < img_shape.1);
    let focal_dis = 0.5 / (fov / 2.0).to_radians().tan();
    let (screen_width, screen_height) = if img_shape.0 > img_shape.1 {
        (img_shape.0 as f32 / img_shape.1 as f32, 1f32)
    } else {
        (1f32, img_shape.1 as f32 / img_shape.0 as f32)
    };
    let x = ((ix as f32 + 0.5 + dx) / img_shape.0 as f32 - 0.5) * screen_width;
    let y = (0.5 - (iy as f32 + 0.5 + dy) / img_shape.1 as f32) * screen_height;
    let z = focal_dis;
    // flip x making right-handed
    let mut dir = [-x, y, z];
    let mut org = [0.0, 0.0, 0.0];
    use del_geo_core::mat4_col_major;
    dir = mat4_col_major::transform_direction(transform_camlcl2world, &dir);
    org = transform_vector(transform_camlcl2world, &org).unwrap();
    (org, dir)
}

pub fn gen_ray_lens(
    (lensx, lensy): (f32, f32),
    (lens_rad, focal_dis): (f32, f32),
    (ix, iy): (usize, usize),
    (dx, dy): (f32, f32),
    img_shape: (usize, usize),
    fov: f32,
    transform_camlcl2world: &[f32; 16],
) -> ([f32; 3], [f32; 3]) {
    assert!(ix < img_shape.0 && iy < img_shape.1);
    let focal_len = 0.5 / (fov / 2.0).to_radians().tan();
    let (screen_width, screen_height) = if img_shape.0 > img_shape.1 {
        (img_shape.0 as f32 / img_shape.1 as f32, 1f32)
    } else {
        (1f32, img_shape.1 as f32 / img_shape.0 as f32)
    };
    let x = ((ix as f32 + 0.5 + dx) / img_shape.0 as f32 - 0.5) * screen_width;
    let y = (0.5 - (iy as f32 + 0.5 + dy) / img_shape.1 as f32) * screen_height;
    let z = focal_len;
    // flip x making right-handed
    let mut dir = [-x, y, z];
    let mut org = [0.0, 0.0, 0.0];

    use del_geo_core::mat4_col_major;
    dir = mat4_col_major::transform_direction(transform_camlcl2world, &dir);
    org = transform_vector(transform_camlcl2world, &org).unwrap();

    if lens_rad > 0. {
        let mut pdisk = sampling::sample_uni_disk_concentric(&[lensx, lensy]);
        pdisk = del_geo_core::vec2::scale(&pdisk, lens_rad);

        let mut plens = vec3::add(&org, &dir);
        plens[0] += pdisk[0];
        plens[1] += pdisk[1];

        let ft = focal_dis;
        let nrm_dir = vec3::normalize(&dir);
        let pfocus = vec3::axpy(ft, &nrm_dir, &org);

        dir = vec3::sub(&pfocus, &plens);
        org = plens;
    }

    (org, dir)
}

/// pos: camera pos
pub fn matrix_w2v(pos: &[Real; 3], view: &[Real; 3]) -> ([Real; 3], [Real; 16]) {
    let (up, view2world) = matrix_v2w(view);
    let mut world2view = del_geo_core::mat4_col_major::transpose(&view2world);
    world2view[0 + 3 * 4] = -pos[0];
    world2view[1 + 3 * 4] = -pos[1];
    world2view[2 + 3 * 4] = -pos[2];

    (up, world2view)
}

pub fn matrix_v2w(view: &[Real; 3]) -> ([Real; 3], [Real; 16]) {
    let forward = &mut vec3::normalize(view);
    let mut up = [0., 1., 0.];
    let equal = vec3::dot(forward, &up) < 1e-7;
    if equal {
        up[2] = 0.1;
        up = vec3::normalize(&up);
    }

    let right = vec3::cross(&up, forward);
    let right = vec3::normalize(&right);
    up = vec3::cross(forward, &right);
    up = vec3::normalize(&up);

    // formatter will not work if use matrix-style initializing
    (
        up,
        [
            right[0], right[1], right[2], 0., up[0], up[1], up[2], 0., forward[0], forward[1],
            forward[2], 0., 0., 0., 0., 1.,
        ],
    )
}

pub fn project_matrix(n: Real, f: Real, fov: Real) -> Matrix4<Real> {
    #[rustfmt::skip]
    let p = Matrix4::new(
        1., 0., 0., 0., 
        0., 1., 0., 0., 
        0., 0., f/(f-n), -f*n/(f-n), 
        0., 0., 1., 0.,
    );

    let invfov = 1. / (fov.to_radians() / 2.).tan();
    #[rustfmt::skip]
    let s = Matrix4::new(
        invfov, 0., 0., 0., 
        0., invfov, 0., 0., 
        0., 0., 1.,0. , 
        0., 0., 0., 1.,
    );

    s * p
}

#[test]
fn test_gen_ray_lens() {
    use crate::cam;

    let lens_rad = 0.05;
    let focal_dis = 2.0;
    let fov = 60.0;

    let w = 5;
    let h = 5;
    let iw = 2;
    let ih = 2;

    let campos = [0.0, 0., 2.86];
    let view = [0., 0., -1.];
    let mut v2w = cam::matrix_v2w(&view).1;
    // concat translation
    v2w[0 + 3 * 4] = campos[0];
    v2w[1 + 3 * 4] = campos[1];
    v2w[2 + 3 * 4] = campos[2];

    let lensx: Real = 0.;
    let lensy: Real = 0.;

    let (ray_org, ray_dir) = cam::gen_ray_lens(
        (lensx, lensy),
        (lens_rad, focal_dis),
        (iw, ih),
        (0., 0.),
        (w, h),
        fov,
        &v2w,
    );

    dbg!(ray_org);
    dbg!(ray_dir);
    let pfocus = vec3::add(&ray_org, &ray_dir);
    dbg!(pfocus);
}
