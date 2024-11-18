use rs_sampler::{
    self, haltonsampler::HaltonSampler, sampler::RandomStrategy, sampler::Sampler,
    sobolsampler::SobolSampler,
};
use std::collections::HashMap;
use std::fs::File;
use std::io::{Error, Write};

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
    let indices = 64;

    // 1st dimension
    for c in 0..indices {
        s.set_i(c);
        s.set_dim(0);
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
        s.restore();
    }

    s.restore();

    // higher dimension
    for c in 0..indices {
        s.set_i(c);
        s.set_dim(63);
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
        s.restore();
    }

    (cols, indices)
}

#[test]
fn test_one_dim_sobol_sampler() {
    let mut s = SobolSampler::new();
    let mut enteries: HashMap<usize, CSVEntery> = HashMap::new();
    let (cols, _indeces) = gen_samples(&mut s, &mut enteries);
    let _ = export_csv(
        "target/sobol_test.csv",
        enteries.len() / cols,
        cols,
        &enteries,
    );
}

#[test]
fn test_halton_sampler() {
    let mut s = HaltonSampler::new_randomized(RandomStrategy::PermuteDigits);
    let mut enteries: HashMap<usize, CSVEntery> = HashMap::new();
    let (cols, _indeces) = gen_samples(&mut s, &mut enteries);
    let _ = export_csv(
        "target/halton_test.csv",
        enteries.len() / cols,
        cols,
        &enteries,
    );
}
