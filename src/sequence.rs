use std::fs::File;
use std::io::prelude::*;
use std::io::{self, BufReader};
use crate::scoring_matrix::ScoringMatrix;

pub fn calc_score(
    _alignment: &Alignment,
    top: &Sequence,
    bot: &Sequence,
    score_mat: &ScoringMatrix,
    gip: i32
) -> i32 {
    let mut t = 0;
    let mut b = 0;
    let mut score = 0;
    for dir in _alignment.directions.iter() {
        match dir {
            DirLinear::Up => {
                t += 1;
                score += gip;
            },
            DirLinear::Right => {
                b += 1;
                score += gip;
            },
            DirLinear::Diag => {
                score += score_mat.m[bot.ix[b]][top.ix[t]];
                t += 1;
                b += 1;
            },
            DirLinear::None => panic!("Invalid direction"),
        };
    }
    //println!("\n --> score: {} # t:{} # b:{} # top:{} # bot:{}", score, t, b, top.len(), bot.len());
    return score;
}

pub struct Sequence {
    pub name: String,
    pub ix: Vec<usize>,
}

pub fn empty_sequence() -> Sequence {
    Sequence {
        name: "".to_string(),
        ix: Vec::new(),
    }
}

impl Sequence {
    pub fn len(&self) -> usize {
        self.ix.len()
    }

    pub fn first_half(&self) -> Sequence {
        let mut name = self.name.to_owned();

        let h = self.ix.len() / 2;

        name.push_str(&format!("[-{}]", h));

        let ix = self.ix[..h].to_vec();

        Sequence { name, ix }
    }

    pub fn second_half(&self) -> Sequence {
        let mut name = self.name.to_owned();

        let h = self.ix.len() / 2;
        name.push_str(&format!("[{}-]", h));
        let ix = self.ix[h..].to_vec();

        Sequence { name, ix }
    }

    pub fn reverse(&self) -> Sequence {
        let mut name = self.name.to_owned();
        name.push_str("[Reversed]");

        let mut ix = self.ix.to_owned();
        ix.reverse();

        Sequence { name, ix }
    }

    pub fn add_sequence(&mut self, seq: &Sequence) {
        self.name.push_str("-");
        self.name.push_str(&seq.name);
        for &elem in &seq.ix {
            self.ix.push(elem);
        }
    }

    pub fn get_subsequence(&self, b: usize, e: usize) -> Sequence {
        let mut name = self.name.to_owned();
        name.push_str(&format!("[{}-{}]", b, e));
        let ix = if e < b {
            Vec::new()
        } else {
            self.ix[b..e].to_vec()
        };

        Sequence { name, ix }
    }
}

pub fn read(filename: &str, c2i: &[i32; 256]) -> io::Result<Sequence> {
    let mut f = BufReader::new(File::open(filename)?);

    let mut line = String::new();
    f.read_line(&mut line)?;
    let first_line: Vec<&str> = line.split("|").collect();
    let name = first_line[1].to_string();

    let mut ix: Vec<usize> = Vec::new();
    for line in f.lines() {
        let sl: String = line.unwrap();
        for elem in sl.chars() {
            ix.push(c2i[elem as usize] as usize);
        }
    }

    Ok(Sequence { name, ix })
}

// ------------------------------------------------------------------------------------------------

#[derive(Copy, Clone, PartialEq)]
#[repr(u8)]
pub enum DirLinear {
    None,
    Up,
    Right,
    Diag,
}

pub fn dirlineal2char(dir: &DirLinear) -> char {
    return match dir {
        DirLinear::None => '-',
        DirLinear::Up => 'U',
        DirLinear::Right => 'R',
        DirLinear::Diag => 'D'
    }
}

pub struct Alignment {
    pub top: Sequence,
    pub bot: Sequence,
    pub directions: Vec<DirLinear>,
}

pub fn empty_alignment() -> Alignment {
    Alignment {
        top: empty_sequence(),
        bot: empty_sequence(),
        directions: Vec::new(),
    }
}

impl Alignment {
    pub fn add_alignment(&mut self, a: &Alignment) {
        self.top.add_sequence(&a.top);
        self.bot.add_sequence(&a.bot);
        for dir in &a.directions {
            self.directions.push(*dir);
        }
    }

    pub fn len(&self) -> usize {
        return self.directions.len();
    }

    pub fn as_string(&self) -> String {
        return self.directions.iter().map(|d| dirlineal2char(d)).collect();
    }
}
