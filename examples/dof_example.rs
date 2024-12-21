use rs_sampler::{envmap::pixel2texpair, haltonsampler::HaltonSampler, sampler::Sampler};

type Real = f32;
type Rgb = image::Rgb<Real>;

fn gen_spheres() -> Vec<([Real; 3], Real)> {
    let mut vec: Vec<([Real; 3], Real)> = Vec::new();

    let nsphere = 3;

    for i in 0..nsphere {
        let c = [i as Real, i as Real, 0.];
        let r = 0.5;
        vec.push((c, r));
    }

    vec
}

fn main() {
    use del_geo_core::vec3;
    use image::Pixel;
    use itertools::Itertools;
    use rayon::iter::IndexedParallelIterator;
    use rayon::iter::IntoParallelRefMutIterator;
    use rayon::iter::ParallelIterator;

    let w = 512;
    let h = 512;
    let nsamples = 64;

    let len_rad = 0.0125;
    let focal_dis = 0.05;
    let fov = 20.0;

    let transform_env = [
        -0.386527, 0., 0.922278, 0., -0.922278, 0., -0.386527, 0., 0., 1., 0., 0., 0., 0., 0., 1.,
    ];
    let transform_env: [f32; 16] = {
        let m = nalgebra::Matrix4::<f32>::from_column_slice(&transform_env);
        let m = m.try_inverse().unwrap();
        // let transform_env = del_geo_core::mat4_col_major::try_inverse(&transform_env).unwrap();
        m.as_slice().try_into().unwrap()
    };
    let transform_cam_lcl2glbl = del_geo_core::mat4_col_major::from_translate(&[0., 0., -5.]);
    let spheres = gen_spheres();

    let mut img = vec![*Rgb::from_slice(&[0.; 3]); w * h];
    let iter = |i_pix: usize, pix: &mut Rgb| {
        let iw = i_pix % w;
        // top-left to bottom-left
        let ih = w - 1 - i_pix / w;
        let tex = pixel2texpair(iw, ih, w, h);

        let (ray_org, ray_dir) = del_raycast_core::cam_pbrt::cast_ray(
            (iw, ih),
            (0., 0.),
            (w, h),
            fov,
            transform_cam_lcl2glbl,
        );

        let mut result = [0.; 3];
        let mut neart = Real::INFINITY;
        let mut cntr = [0.; 3];
        for sphere in &spheres {
            let hit_res = del_geo_nalgebra::sphere::intersection_ray(
                &nalgebra::Vector3::<f32>::from(sphere.0),
                0.7,
                &nalgebra::Vector3::<f32>::from(ray_org),
                &nalgebra::Vector3::<f32>::from(ray_dir),
            );

            if let Some(t) = hit_res {
                if t < neart {
                    neart = t;
                    cntr = sphere.0;
                }
            }
        }
        let hit_pos = vec3::axpy::<f32>(neart, &ray_dir, &ray_org);
        let hit_nrm = vec3::normalize(&vec3::sub(&hit_pos, &cntr));

        let env =
            del_geo_core::mat4_col_major::transform_homogeneous(&transform_env, &hit_nrm).unwrap();

        result = env;

        // for i in 0..nsamples {
        //     todo!()
        // }

        pix.0 = result;
    };

    img.par_iter_mut()
        .enumerate()
        .for_each(|(i_pix, pix)| iter(i_pix, pix));

    use image::codecs::hdr::HdrEncoder;
    let file_ms = std::fs::File::create("target/dof.hdr").unwrap();
    let enc = HdrEncoder::new(file_ms);
    let _ = enc.encode(&img, w, h);
}
