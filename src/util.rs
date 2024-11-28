use std::collections::HashMap;
use std::fs::File;
use std::io::{Error, Write};

#[derive(Default)]
pub struct CSVEntery {
   pub val: String,
}

impl Clone for CSVEntery {
    fn clone(&self) -> Self {
        CSVEntery::default()
    }
}

pub fn export_csv(
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

