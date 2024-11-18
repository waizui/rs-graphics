use std::fs::File;
use std::io::{Error, Write};
use std::vec;
use rs_sampler::{self, sampler::Sampler, sobolsampler::SobolSampler};

#[derive(Default)]
struct CSVEntery {
    r: usize,
    c: usize,
    val: usize,
}

impl Clone for CSVEntery {
    fn clone(&self) -> Self {
        CSVEntery::default()
    }
}

fn export_csv(path: &str, mut a: &[CSVEntery]) -> Result<(), Error> {
    let mut file = File::create(path)?;

    let rows = (a.len() as f32).sqrt() as usize;

    for r in 0..rows {
        for i in 0..rows {
            let e = &a[r * rows + i];
            if e.val != 0 {
                let _ = file.write(format!("{},", e.val).as_bytes());
            } else {
                let _ = file.write(b"0,");
            }
        }

        let _ = file.write(b"\n");
    }

    file.flush()?;
    Ok(())
}

#[test]
fn test_sobolsampler() {
    let mut s = SobolSampler::new();
    let rows = 16;
    let c_i = rows * rows;

    let mut enteries = vec![CSVEntery::default(); c_i];
    // let rand = StdRng::seed_from_u64(1u64).gen_range(0..c_i);
    for i in 0..c_i {
        s.index = i;
        let v: [f32; 2] = s.get2d();

        let r: usize = (v[0] * (rows as f32)).floor() as usize;
        let c: usize = (v[1] * (rows as f32)).floor() as usize;
        let cell_i = r * rows + c;
        enteries[cell_i].r = r;
        enteries[cell_i].c = c;
        enteries[cell_i].val += 1;
        s.restore();
    }

    let _ = export_csv("target/one_dim_sobol_test.csv", &enteries[0..]);

    let mut s = SobolSampler::new();
    let mut enteries = vec![CSVEntery::default(); c_i];
    // let rand = StdRng::seed_from_u64(1u64).gen_range(0..c_i);
    for i in 0..c_i {
        s.index = i;
        s.dim = 64;
        let v: [f32; 2] = s.get2d();

        let r: usize = (v[0] * (rows as f32)).floor() as usize;
        let c: usize = (v[1] * (rows as f32)).floor() as usize;
        let cell_i = r * rows + c;
        enteries[cell_i].r = r;
        enteries[cell_i].c = c;
        enteries[cell_i].val += 1;
        s.restore();
        // let entery = CSVEntery {};
        // enteries.push(entery);
    }

    let _ = export_csv("target/mul_dim_sobol_test.csv", &enteries[0..]);
}
