use crate::baseline::baseline_get_score;
use crate::quadratic_space::qs;
use crate::scoring_matrix::ScoringMatrix;
use crate::sequence::{Alignment, Sequence, calc_score};

fn combine_scores(col_d: &Vec<i32>, col_r: &Vec<i32>) -> (i32, usize) {
    //Best score on central column
    let mut best_center_score = col_d[0] + col_r[col_d.len() - 1];
    let mut index_best_score = 0;
    for i in 1..col_d.len() {
        let ns = col_d[i] + col_r[col_d.len() - 1 - i];
        if ns > best_center_score {
            best_center_score = ns;
            index_best_score = i;
        }
    }
    (best_center_score, index_best_score)
}

pub fn hirschberg(
    top: &Sequence,
    bot: &Sequence,
    score_mat: &ScoringMatrix,
    gip: i32,
    threshold: usize,
) -> (i32, Alignment) {
    let left = bot.first_half();
    let right = bot.second_half();
    let top_reversed = top.reverse();
    let right_reversed = right.reverse();

    let left_col = baseline_get_score(&top, &left, &score_mat, gip);
    let right_col = baseline_get_score(&top_reversed, &right_reversed, &score_mat, gip);
    let (best_score, index_best_score) = combine_scores(&left_col, &right_col);

    let top1 = top.get_subsequence(0, index_best_score);
    let top2 = top.get_subsequence(index_best_score, top.len());

    let (sco_l, mut d) = if left.len() <= 3 || left.len() * top1.len() < threshold {
        qs(&top1, &left, &score_mat, gip)
    } else {
        hirschberg(&top1, &left, &score_mat, gip, threshold)
    };
    debug_assert!(calc_score(&d, &top1, &left, &score_mat, gip) == sco_l);

    let (sco_r, r) = if right.len() <= 3 || right.len() * top2.len() < threshold {
        qs(&top2, &right, &score_mat, gip)
    } else {
        hirschberg(&top2, &right, &score_mat, gip, threshold)
    };
    debug_assert!(calc_score(&r, &top2, &right, &score_mat, gip) == sco_r);

    debug_assert!(best_score == sco_l+sco_r);

    d.add_alignment(&r);
    debug_assert!(calc_score(&d, &top, &bot, &score_mat, gip) == sco_l+sco_r);

    (best_score, d)
}
