use crate::scoring_matrix::ScoringMatrix;
use crate::sequence::Sequence;

pub fn baseline_get_score(top: &Sequence, bot: &Sequence, score_mat: &ScoringMatrix, gip: i32) -> Vec<i32> {
    let mut scores: Vec<i32> = vec![0; top.len() + 1];

    scores[0] = 0;

    for t in 1..scores.len() {
        scores[t] = scores[t - 1] + gip;
    }

    for b in 1..bot.len() + 1 {
        let am = score_mat.m[bot.ix[b - 1]];

        let mut pre_old_score = scores[0];
        scores[0] = scores[0] + gip;
        let mut pre_score = scores[0];
        for t in 1..scores.len() {
            let score_diag = pre_old_score + am[top.ix[t - 1]];
            pre_old_score = scores[t];
            let score_left = pre_old_score + gip;
            let score_up = pre_score + gip;

            pre_score = if score_diag > score_left {
                if score_diag > score_up {
                    score_diag
                } else {
                    score_up
                }
            } else {
                if score_left > score_up {
                    score_left
                } else {
                    score_up
                }
            };
            scores[t] = pre_score;
        }
    }

    scores
}

