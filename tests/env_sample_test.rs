use del_raycast_core::io_pfm::PFM;

type Image = PFM;

fn calc_gray_scale(path: &str) -> Image {
    let pfm = PFM::read_from(path).unwrap();

    if pfm.channels < 3 {
        panic!("pfm channel invalid");
    }

    let mut gray_scale = PFM {
        data: vec![0.; pfm.w * pfm.h],
        w: pfm.w,
        h: pfm.h,
        channels: 1,
        little_endian: pfm.little_endian,
    };

    let channels = pfm.channels;
    for h in 0..pfm.h {
        for w in 0..pfm.w {
            let p_i = channels * (h * pfm.w + w);
            let r = pfm.data[p_i];
            let g = pfm.data[p_i + 1];
            let b = pfm.data[p_i + 2];
            let gray = (0.2126 * r + 0.7152 * g + 0.0722 * b).clamp(0., 1.);
            gray_scale.data[h * pfm.w + w] = gray;
        }
    }

    gray_scale
}

fn calc_cdf_inv(img: &Image) -> [Image; 2] {
    todo!()
}

#[test]
fn test_env_sample() {
    let img = calc_gray_scale("asset/envmap.pfm");
    // let _ = img.save("target/envmap_gray.pfm");
    // let cdf_inv = calc_cdf_inv(&img);
}
