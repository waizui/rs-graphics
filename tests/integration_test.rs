use rs_sampler::util::{export_csv, CSVEntery};

use rs_sampler::{
    self, haltonsampler::HaltonSampler, sampler::RandomStrategy, sampler::Sampler,
    sobolsampler::SobolSampler,
};
use std::collections::HashMap;
use std::time::Instant;

fn gen_samples<T>(s: &mut T, map: &mut HashMap<usize, CSVEntery>) -> (usize, usize)
where
    T: Sampler<f32>,
{
    let cols = 4;
    let indices = 256;

    // 1st dimension
    for c in 0..indices {
        s.set_dim(0);
        s.set_i(c);
        let v: [f32; 2] = s.get2d();

        map.insert(
            c * cols,
            CSVEntery {
                val: v[0].to_string(),
            },
        );
        map.insert(
            c * cols + 1,
            CSVEntery {
                val: v[1].to_string(),
            },
        );
    }

    s.restore();

    //higher dimension
    for c in 0..indices {
        s.set_dim(63);
        s.set_i(c);
        let v: [f32; 2] = s.get2d();

        map.insert(
            c * cols + 2,
            CSVEntery {
                val: v[0].to_string(),
            },
        );
        map.insert(
            c * cols + 3,
            CSVEntery {
                val: v[1].to_string(),
            },
        );
    }

    (cols, indices)
}

fn perf_samples<T>(s: &mut T, loops: usize, count: usize, id: &str)
where
    T: Sampler<f32>,
{
    let start = Instant::now();

    for _ in 0..loops {
        for c in 0..count {
            s.set_i(c);
            let v: [f32; 2] = s.get2d();
            let _ = drop(v);
        }
        s.restore();
    }

    let elapsed = start.elapsed().as_millis();
    println!("{} perf:{} ms", id, elapsed);
}

#[test]
fn test_sobol_sampler() {
    let mut s = SobolSampler::new_randomized(RandomStrategy::None);
    let mut enteries: HashMap<usize, CSVEntery> = HashMap::new();
    let (cols, _indeces) = gen_samples(&mut s, &mut enteries);
    let _ = export_csv(
        "target/sobol_test_norandom.csv",
        enteries.len() / cols,
        cols,
        &enteries,
    );

    let mut s = SobolSampler::new_randomized(RandomStrategy::PermuteDigits);
    enteries.clear();
    let (cols, _indeces) = gen_samples(&mut s, &mut enteries);
    let _ = export_csv(
        "target/sobol_test_permute.csv",
        enteries.len() / cols,
        cols,
        &enteries,
    );
}

#[test]
fn test_halton_sampler() {
    let mut s = HaltonSampler::new_randomized(RandomStrategy::None);
    let mut enteries: HashMap<usize, CSVEntery> = HashMap::new();
    let (cols, _indeces) = gen_samples(&mut s, &mut enteries);
    let _ = export_csv(
        "target/halton_test_norandom.csv",
        enteries.len() / cols,
        cols,
        &enteries,
    );

    let mut s = HaltonSampler::new_randomized(RandomStrategy::PermuteDigits);
    enteries.clear();
    let (cols, _indeces) = gen_samples(&mut s, &mut enteries);
    let _ = export_csv(
        "target/halton_test_permute.csv",
        enteries.len() / cols,
        cols,
        &enteries,
    );
}

#[test]
fn test_halton_perf() {
    let mut s = HaltonSampler::new_randomized(RandomStrategy::PermuteDigits);
    perf_samples(&mut s, 128, 256, "halton");
}

#[test]
fn test_sobol_perf() {
    let mut s = SobolSampler::new_randomized(RandomStrategy::None);
    perf_samples(&mut s, 128, 256, "sobol");
}
