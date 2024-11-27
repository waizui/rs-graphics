use rs_sampler::{
    self, haltonsampler::HaltonSampler, sampler::RandomStrategy, sampler::Sampler,
    sobolsampler::SobolSampler,
};
use std::collections::HashMap;
use std::fs::File;
use std::io::{Error, Write};
use std::time::Instant;

#[derive(Default)]
struct CSVEntery {
    val: String,
}

impl Clone for CSVEntery {
    fn clone(&self) -> Self {
        CSVEntery::default()
    }
}

fn export_csv(
    path: &str,
    rows: usize,
    cols: usize,
    a: &HashMap<usize, CSVEntery>,
) -> Result<(), Error> {
    let mut file = File::create(path)?;
    for rol in 0..rows {
        for col in 0..cols {
            let cell_i = rol * cols + col;
            if let Some(e) = a.get(&cell_i) {
                let _ = file.write(format!("{},", e.val).as_bytes());
            } else {
                let _ = file.write(b",");
            }
        }
        let _ = file.write(b"\n");
    }

    file.flush()?;
    Ok(())
}

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
