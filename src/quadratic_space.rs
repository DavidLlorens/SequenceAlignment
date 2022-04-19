use crate::scoring_matrix::ScoringMatrix;
use crate::sequence::{Alignment, DirLinear, Sequence};

pub fn qs(top: &Sequence, bot: &Sequence, score_mat: &ScoringMatrix, gip: i32) -> (i32, Alignment) {
    let mut old_score: Vec<i32> = vec![0; top.len() + 1];
    let mut new_score: Vec<i32> = vec![0; top.len() + 1];
    let mut bp: Vec<Vec<DirLinear>> = vec![vec![DirLinear::None; top.len() + 1]; bot.len() + 1];

    bp[0][0] = DirLinear::None;
    old_score[0] = 0;

    for t in 1..old_score.len() {
        bp[0][t] = DirLinear::Up;
        old_score[t] = old_score[t - 1] + gip;
    }

    for b in 1..bot.len() + 1 {
        let am = score_mat.m[bot.ix[b - 1]];

        new_score[0] = old_score[0] + gip;
        let mut pre_old_score = old_score[0];
        let mut pre_score = new_score[0];
        let bpb = &mut bp[b];
        bpb[0] = DirLinear::Right;
        for t in 1..old_score.len() {
            let score_diag = pre_old_score + am[top.ix[t - 1]];
            pre_old_score = old_score[t];
            let score_left = pre_old_score + gip;
            let score_up = pre_score + gip;

            let direction: DirLinear;
            pre_score = if score_diag > score_left {
                if score_diag > score_up {
                    direction = DirLinear::Diag;
                    score_diag
                } else {
                    direction = DirLinear::Up;
                    score_up
                }
            } else {
                if score_left > score_up {
                    direction = DirLinear::Right;
                    score_left
                } else {
                    direction = DirLinear::Up;
                    score_up
                }
            };
            bpb[t] = direction;
            new_score[t] = pre_score;
        }
        let tmp = old_score;
        old_score = new_score;
        new_score = tmp;
    }

    let b_e = bot.len();
    let mut b_b = b_e;

    let t_e = top.len();
    let mut t_b = t_e;

    let mut dirs: Vec<DirLinear> = Vec::with_capacity(b_e + t_e);
    while bp[b_b][t_b] != DirLinear::None {
        dirs.push(bp[b_b][t_b]);
        match bp[b_b][t_b] {
            DirLinear::Diag => {
                b_b -= 1;
                t_b -= 1;
            }
            DirLinear::Up => t_b -= 1,
            DirLinear::Right => b_b -= 1,
            _ => {}
        }
    }
    dirs.reverse();
    debug_assert!(b_b==0 && t_b==0, "b_b=={} # t_b=={}", b_b, t_b);
    //println!("---> {} {} {}", dirs.len(), top.len(), bot.len());
    (
        old_score[top.len()],
        Alignment {
            top: top.get_subsequence(t_b, t_e),
            bot: bot.get_subsequence(b_b, b_e),
            directions: dirs,
        },
    )
}
