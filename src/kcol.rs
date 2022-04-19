use crate::quadratic_space::qs;
use crate::scoring_matrix::ScoringMatrix;
use crate::sequence::{Alignment, Sequence, empty_alignment, calc_score};

fn get_intercol_segments(
    cols: &Vec<usize>,
    stored_cols: &Vec<Vec<usize>>,
    top: &Sequence,
    bot: &Sequence,
) -> Vec<(Sequence, Sequence)> {
    let mut res: Vec<(Sequence, Sequence)> = Vec::new();
    let mut t_e = top.len();
    //Debug.Assert(_cols[posCols] == posStoredCols);
    for pos_cols in (1..cols.len()).rev() {
        let mut t_b = stored_cols[pos_cols][t_e];
        if pos_cols == 1 {
            t_b = 0;
        }
        //println!("t({},{}) b({},{})", t_b, t_e, cols[pos_cols - 1], cols[pos_cols]);
        res.push((
            top.get_subsequence(t_b, t_e),
            bot.get_subsequence(cols[pos_cols - 1], cols[pos_cols]),
        ));
        t_e = t_b;
    }
    //println!("Score: {}", stored_cols[stored_cols.len() - 1][top.len()]);
    res.reverse();
    return res;
}

pub fn kcol(
    top: &Sequence,
    bot: &Sequence,
    score_mat: &ScoringMatrix,
    gip: i32,
    num_regions: usize,
    threshold: usize,
) -> (i32, Alignment) {
    let num_stored_cols: usize = num_regions + 1;
    let mut cols: Vec<usize> = vec![0; num_stored_cols];
    cols[0] = 0;
    cols[num_stored_cols - 1] = bot.len();
    let s: f64 = bot.len() as f64 / (num_regions as f64);
    for i in 1..num_stored_cols - 1 {
        cols[i] = (s * i as f64) as usize;
    }
    let mut stored_cols: Vec<Vec<usize>> = vec![vec![0; top.len() + 1]; num_stored_cols]; //new int[num_stored_cols][];
    let mut scores: Vec<(i32, usize)> = vec![(0, 0); top.len() + 1];

    scores[0] = (0, 0);

    for t in 1..scores.len() {
        scores[t] = (scores[t - 1].0 + gip, t);
    }

    let mut pos_in_cols = 1;

    for b in 1..bot.len() + 1 {
        let am = score_mat.m[bot.ix[b - 1]];

        let mut pre_old_score = scores[0];
        scores[0] = (scores[0].0 + gip, scores[0].1);
        let mut pre_score = scores[0];
        for t in 1..scores.len() {
            let score_diag = (pre_old_score.0 + am[top.ix[t - 1]], pre_old_score.1);
            pre_old_score = scores[t];
            let score_left = pre_old_score.0 + gip;
            let score_up = pre_score.0 + gip;

            pre_score = if score_diag.0 > score_left {
                if score_diag.0 > score_up {
                    score_diag
                } else {
                    (score_up, pre_score.1)
                }
            } else {
                if score_left > score_up {
                    (score_left, pre_old_score.1)
                } else {
                    (score_up, pre_score.1)
                }
            };
            scores[t] = pre_score;
        }
        if b == cols[pos_in_cols] {
            let scp = &mut stored_cols[pos_in_cols];
            for i in 0..scores.len() {
                scp[i] = scores[i].1;
                scores[i].1 = i;
            }
            pos_in_cols += 1;
        }
    }

    let mut align = empty_alignment();

    for (top2, bot2) in get_intercol_segments(&cols, &stored_cols, top, bot) {
        let (_score, al) =
            if bot2.len() <= num_regions || bot2.len() * top2.len() < threshold {
                qs(&top2, &bot2, &score_mat, gip)
            } else {
                kcol(&top2, &bot2, &score_mat, gip, num_regions, threshold)
            };
        debug_assert!(calc_score(&al, &top2, &bot2, &score_mat, gip) == _score);
        align.add_alignment(&al);
    }
    debug_assert!(calc_score(&align, &top, &bot, &score_mat, gip) == scores[top.len()].0);
    (scores[top.len()].0, align)
}
