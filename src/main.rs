extern crate getopts;
use std::fs::File;
use std::io::Write;

use getopts::Options;

use std::env;
use std::time::Instant;

mod scoring_matrix;

mod sequence;

mod baseline;
use baseline::baseline_get_score;

mod quadratic_space;
use quadratic_space::qs;

mod hirschberg;
use hirschberg::hirschberg;

mod kcol;
use kcol::kcol;

const VERSION: &str = env!("CARGO_PKG_VERSION");

fn print_usage(program: &str, opts: Options) {
    println!("Program version: {}\n", VERSION);
    let brief = format!("Usage {} [program options] MODE [mode options]", program);
    print!("{}", opts.usage(&brief));
}

const DEFAULT_TOP: &str = "TestData/A2ASS6.fasta";
const DEFAULT_BOT: &str = "TestData/Q8WZ42.fasta";
const DEFAULT_SCORES: &str = "TestData/blosum62.txt";
const DEFAULT_GIP: i32 = -10;

const DEFAULT_THRESHOLD: usize = 100000;
const DEFAULT_NUM_REGIONS: usize = 32;

fn main() {
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();

    let mut opts = Options::new();
    opts.optflag("h", "help", "This help");
    opts.optopt("b", "bottom", &format!("file with the bottom sequence (default: {})", DEFAULT_BOT), "FILE");
    opts.optopt("t", "top", &format!("file with the top sequence (default: {})", DEFAULT_TOP), "FILE");
    opts.optopt("s", "scores", &format!("The file with the score matrix (default: {})", DEFAULT_SCORES), "FILE");
    opts.optopt("g", "gip", &format!("The gap insertion penalty (default: {})", DEFAULT_GIP), "INT");
    opts.optopt("a", "align", &"file where the alignment will be saved (default: not saved)", "FILE");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => { m }
        Err(f) => { panic!("{}", f.to_string()) }
    };

    if matches.opt_present("h") {
        print_usage(&program, opts);
        return ;
    }

    let mode = if !matches.free.is_empty() {
        matches.free[0].clone()
    } else {
        print_usage(&program, opts);
        return;
    };

    let gip = matches.opt_get_default("g", DEFAULT_GIP).unwrap();
    let score_file : &str = &matches.opt_get_default("s", DEFAULT_SCORES.to_string()).unwrap();
    let score_mat = scoring_matrix::read(score_file).unwrap();
    let top_file : &str = &matches.opt_get_default("t", DEFAULT_TOP.to_string()).unwrap();
    let top = sequence::read(top_file, &score_mat.c2i).unwrap();
    let bot_file : &str = &matches.opt_get_default("b", DEFAULT_BOT.to_string()).unwrap();
    let bot = sequence::read(bot_file, &score_mat.c2i).unwrap();
    let align_file : &str = &matches.opt_get_default("a", "".to_string()).unwrap();

    let mut alignment = sequence::empty_alignment();

    match mode.as_str() {
        "info" => {
            println!("Program version: {}\n", VERSION);
            println!("Symbols in the score matrix: {}", score_mat.symbols);
            println!("\nTop sequence:\n    - Name: {}\n    - Length: {}", top.name, top.ix.len());
            println!("Bottom sequence:\n    - Name: {}\n    - Length: {}", bot.name, bot.ix.len());
        }
        "bl" => {
            let instant = Instant::now();
            let score = baseline_get_score(&top, &bot, &score_mat, gip);
            println!("\nBaseline score: {}", score[top.len()]);
            println!("Time: {}", instant.elapsed().as_secs_f64());
        }
        "qs" => {
            let instant = Instant::now();
            let (score, _alignment) = qs(&top, &bot, &score_mat, gip);
            alignment = _alignment;
            println!("\nQuadratic-Alg score: {}", score);
            println!("Quadratic-Alg alignment length: {}", alignment.len());
            println!("Time: {}", instant.elapsed().as_secs_f64());
            debug_assert!(sequence::calc_score(&alignment, &top, &bot, &score_mat, gip) == score);
        }
        "hb" => {
            let threshold = if matches.free.len() > 1 {
                matches.free[1].parse::<usize>().unwrap()
            } else { DEFAULT_THRESHOLD };

            let instant = Instant::now();
            let (score, _alignment) = hirschberg(&top, &bot, &score_mat, gip, threshold);
            alignment = _alignment;
            println!("\nThreshold: {}", threshold);
            println!("Hirschberg score: {}", score);
            println!("Hirschberg alignment length: {}", alignment.len());
            println!("Time: {}", instant.elapsed().as_secs_f64());
            debug_assert!(sequence::calc_score(&alignment, &top, &bot, &score_mat, gip) == score);
        }
        "kcol" => {
            let num_regions = if matches.free.len() > 1 {
                matches.free[1].parse::<usize>().unwrap()
            } else { DEFAULT_NUM_REGIONS };
            let threshold = if matches.free.len() > 2 {
                matches.free[2].parse::<usize>().unwrap()
            } else { DEFAULT_THRESHOLD };

            let instant = Instant::now();
            let (score, _alignment) = kcol(&top, &bot, &score_mat, gip, num_regions, threshold);
            alignment = _alignment;
            println!("Num regions: {}", num_regions);
            println!("Threshold: {}", threshold);
            println!("K-Col score: {}", score);
            println!("K-Col alignment length: {}", alignment.len());
            println!("Time: {}", instant.elapsed().as_secs_f64());
            debug_assert!(sequence::calc_score(&alignment, &top, &bot, &score_mat, gip) == score);
        }
        _ => { panic!("Unknown mode {}", mode); }
    }
    if align_file.len() > 0 {
        let mut output = File::create(align_file).unwrap();
        let _k = output.write(alignment.as_string().as_bytes());
    }
}
