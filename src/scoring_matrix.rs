use std::fs::File;
use std::io::prelude::*;
use std::io::{self, BufReader};

pub struct ScoringMatrix {
    pub c2i: [i32; 256],
    pub symbols: String,
    pub m: [[i32; 24]; 24],
}

pub fn read(filename: &str) -> io::Result<ScoringMatrix> {
    let mut c2i: [i32; 256] = [-1; 256];
    let mut m: [[i32; 24]; 24] = [[0; 24]; 24];

    let mut f = BufReader::new(File::open(filename)?);

    let mut line = String::new();
    f.read_line(&mut line)?;
    let sl: Vec<&str> = line.split_whitespace().collect();
    let symbols: String = sl.join("");

    for (i, c) in symbols.chars().enumerate() {
        c2i[c as usize] = i as i32;
    }

    let mut r = 0;
    for line in f.lines() {
        let sl: String = line.unwrap();
        let kk: Vec<&str> = sl.split_whitespace().collect();
        for (c, elem) in kk.iter().enumerate() {
            m[r][c] = elem.parse::<i32>().unwrap();
        }
        r += 1;
    }

    Ok(ScoringMatrix { c2i, symbols, m })
}
