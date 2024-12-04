use rs_sampler::envmap;
use rs_sampler::{
    envmap::pixel2tex, envmap::tex2pixel, envmap::Image, haltonsampler::HaltonSampler,
    sampler::Sampler,
};

use rs_sampler::pfm::PFM;

type Real = f32;
const PI: Real = std::f32::consts::PI;

fn binary_find_inv_cdf() -> Real {
    todo!()
}

// fn sample_light(grayscale: &Image, envmap: &Image) -> Image {
//     let w = 64;
//     let h = 64;
//
//     let debug = false;
//
//     let mut img = create_image(3, w, h);
//
//     let mut halton = HaltonSampler::new();
//     let itgr = calc_integral_over_grayscale(grayscale);
//
//     let nsamples = 128;
//     for i_w in 0..w {
//         for i_h in 0..h {
//             // screen coord to world
//             let tex = pixel2tex(i_w, i_h, w, h);
//             let nrm = envmap2unitsphere(&tex);
//
//             let mut resut = [0.; 3];
//             for i in 0..nsamples {
//                 Sampler::<Real>::set_i(&mut halton, i + 1);
//                 Sampler::<Real>::set_dim(&mut halton, 0);
//                 let random: [Real; 2] = halton.get2d();
//                 // sampling position
//                 let coord_pfd = calc_sample_coord(random[0], random[1], grayscale, itgr);
//                 let st = coord_pfd.0;
//                 let pdf = coord_pfd.1;
//                 let s_tex_coord = pixel2tex(st[0], st[1], grayscale.w, grayscale.h);
//                 // sampling ray direction
//                 let raydir = envmap2unitsphere(&s_tex_coord);
//                 let radiance = get_color(s_tex_coord[0], s_tex_coord[1], envmap);
//                 if debug {
//                     set_color(&[1., 0., 0.], s_tex_coord[0], s_tex_coord[1], &mut img);
//                 } else {
//                     let costheta = del_geo_core::vec3::dot(&raydir, &nrm);
//                     if costheta > 0. && pdf > 0. {
//                         let sintheta = (1. - raydir[2] * raydir[2]).sqrt();
//                         let mut color = [PI * 2.; 3]; //Jacobian TODO:replace with correct term
//                         color[0] *= radiance[0];
//                         color[1] *= radiance[1];
//                         color[2] *= radiance[2];
//                         del_geo_core::vec3::scale(&mut color, costheta);
//                         del_geo_core::vec3::scale(&mut color, sintheta);
//                         del_geo_core::vec3::scale(&mut color, 1. / pdf);
//                         resut = del_geo_core::vec3::add(&resut, &radiance);
//                     }
//                 }
//             }
//
//             if debug {
//                 let color = get_color(tex[0], tex[1], &img);
//                 if color[0] < 1. {
//                     let tex_color = get_color(tex[0], tex[1], envmap);
//                     set_color(&tex_color, tex[0], tex[1], &mut img);
//                 }
//             } else {
//                 del_geo_core::vec3::scale(&mut resut, 1. / (nsamples as Real));
//                 set_color(&resut, tex[0], tex[1], &mut img);
//             }
//         }
//     }
//
//     img
// }

#[test]
fn test_inverse_cdf() {
    let pfm = PFM::read_from("asset/envmap.pfm").unwrap();
    let w = pfm.w();
    let h = pfm.h();
    let invcdfmap = envmap::calc_inverse_cdf_map(&pfm);
    let mut res = PFM::new(w, h);
    for p_h in 0..h {
        for p_w in 0..w {
            let tex = invcdfmap.get_p(&[p_w, p_h]);
            let u = tex2pixel(tex[0], w);
            let v = tex2pixel(tex[1], h);
            res.set_p(&[u, v], &[1.; 3]);
        }
    }

    let _ = res.save_to("target/invcdf.pfm");
}

#[test]
fn test_env_sample() {
    // let pfm = PFM::read_from("asset/envmap.pfm").unwrap();
    // let env = sample_light(&grayscale, &pfm);
    // let _ = env.save("target/envmap_gray.pfm");
}
