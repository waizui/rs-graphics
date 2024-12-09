use rs_sampler::pfm::PFM;
use rs_sampler::{
    envmap::envmap2unitsphere, envmap::pixel2texpair, envmap::tex2pixel,
    haltonsampler::HaltonSampler, sampler::Sampler,
};

type Real = f32;
const PI: Real = std::f32::consts::PI;
type Rgb = image::Rgb<Real>;

fn main() {
    use image::Pixel;
    use itertools::Itertools;
    use rayon::iter::IndexedParallelIterator;
    use rayon::iter::IntoParallelRefMutIterator;
    use rayon::iter::ParallelIterator;

    let pfm = PFM::read_from("asset/envmap.pfm").unwrap();
    let w = pfm.w;
    let h = pfm.h;

    let rgbdata = pfm
        .data
        .chunks(pfm.channels)
        .map(|chunk| *Rgb::from_slice(&[chunk[0], chunk[1], chunk[2]]))
        .collect_vec();

    let grayscale = rs_sampler::envmap::calc_grayscale(&rgbdata, w, h);
    let itgr = rs_sampler::envmap::calc_integral_over_grayscale(&grayscale, w, h);

    let (marginal_map, conditional_map) = rs_sampler::envmap::calc_inverse_cdf_map(&rgbdata, w, h);

    let mut img = vec![*Rgb::from_slice(&[0.; 3]); w * h];

    let iter = |i_pix: usize, pix: &mut Rgb| {
        let mut halton = HaltonSampler::new();
        let nsamples = 128;
        let i_w = i_pix % w;
        let i_h = i_pix / w;

        let tex = pixel2texpair(i_w, i_h, w, h);
        let nrm = envmap2unitsphere(&tex);

        let mut result = [0.; 3];

        for i in 0..nsamples {
            Sampler::<Real>::set_i(&mut halton, i + 1);
            Sampler::<Real>::set_dim(&mut halton, 0);
            let random: [Real; 2] = halton.get2d();

            let r_w = tex2pixel(random[0], w);
            let r_h = tex2pixel(random[1], h);

            let sampley = marginal_map[r_h][0];
            let samplex = conditional_map[r_w + tex2pixel(sampley, h) * w][0];

            let p_w = tex2pixel(samplex, w);
            let p_h = tex2pixel(sampley, h);

            let mut radiance = rgbdata[p_h * w + p_w].0;

            let ray_dir = envmap2unitsphere(&[samplex, sampley]);
            let costheta = del_geo_core::vec3::dot(&nrm, &ray_dir);

            if costheta < 0. {
                continue;
            }

            let sintheta = (1. - costheta * costheta).sqrt();

            let pdf = grayscale[p_h * w + p_w][0] / itgr;

            del_geo_core::vec3::scale(&mut radiance, costheta);
            del_geo_core::vec3::scale(&mut radiance, sintheta);
            del_geo_core::vec3::scale(&mut radiance, 1. / pdf);

            result = del_geo_core::vec3::add(&result, &radiance);
        }
        del_geo_core::vec3::scale(&mut result, (4.) / nsamples as Real);

        pix.0 = result;
    };

    img.par_iter_mut()
        .enumerate()
        .for_each(|(i_pix, pix)| iter(i_pix, pix));

    use image::codecs::hdr::HdrEncoder;
    let file_ms = std::fs::File::create("target/04_env_light_sampling.hdr").unwrap();
    let enc = HdrEncoder::new(file_ms);
    let _ = enc.encode(&img, w, h);
}
