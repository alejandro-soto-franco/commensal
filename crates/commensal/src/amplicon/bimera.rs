use rayon::prelude::*;

use crate::Unique;

/// Settings for bimera detection, named as the reference implementation names
/// them so a run can be matched against it.
#[derive(Clone, Copy, Debug)]
pub struct BimeraOptions {
    /// A parent must be at least this many times as abundant as the query.
    /// The reference implementation defaults to 2.
    pub min_fold_parent_over_abundance: f64,
    /// A parent must have at least this abundance outright.
    pub min_parent_abundance: u32,
    /// How far a parent may be shifted against the query when matching. The
    /// reference implementation defaults to 16, which is what lets it call
    /// chimeras whose parents differ in length or have an indel.
    pub max_shift: isize,
}

impl Default for BimeraOptions {
    fn default() -> Self {
        BimeraOptions {
            min_fold_parent_over_abundance: 2.0,
            min_parent_abundance: 8,
            max_shift: 16,
        }
    }
}

/// The longest prefix and suffix of `query` that `parent` covers, over every
/// offset within `max_shift`.
///
/// Comparing `query[i]` against `parent[i + d]` for a shift `d` is what lets a
/// parent of a different length, or one with an indel, still supply a side
/// of the join. The best prefix and the best suffix are taken independently,
/// since the two sides of a chimera come from different parents anyway.
fn coverage(query: &[u8], parent: &[u8], max_shift: isize) -> (usize, usize) {
    let n = query.len() as isize;
    let m = parent.len() as isize;
    let (mut best_prefix, mut best_suffix) = (0usize, 0usize);

    for d in -max_shift..=max_shift {
        let mut i = 0isize;
        while i < n {
            let j = i + d;
            if j < 0 || j >= m || parent[j as usize] != query[i as usize] {
                break;
            }
            i += 1;
        }
        best_prefix = best_prefix.max(i.max(0) as usize);

        let mut k = 0isize;
        while k < n {
            let i = n - 1 - k;
            let j = i + d;
            if j < 0 || j >= m || parent[j as usize] != query[i as usize] {
                break;
            }
            k += 1;
        }
        best_suffix = best_suffix.max(k.max(0) as usize);
    }
    (best_prefix, best_suffix)
}

/// Whether one sequence reads as the left of some parent and the right of
/// another, both drawn from `parents`.
///
/// The split runs over every interior position. A parent pair counts only when
/// the two parents differ and each covers its own side, up to `max_shift`
/// positions of offset.
#[must_use]
pub fn is_bimera(query: &[u8], parents: &[&Unique], max_shift: isize) -> bool {
    let n = query.len();
    if n < 2 || parents.len() < 2 {
        return false;
    }

    let mut prefix = Vec::with_capacity(parents.len());
    let mut suffix = Vec::with_capacity(parents.len());
    for p in parents {
        if p.sequence == query {
            // The query is not one of its own parents.
            prefix.push(0);
            suffix.push(0);
            continue;
        }
        let (a, b) = coverage(query, &p.sequence, max_shift);
        prefix.push(a);
        suffix.push(b);
    }

    for s in 1..n {
        let lefts: Vec<usize> = (0..parents.len()).filter(|&k| prefix[k] >= s).collect();
        if lefts.is_empty() {
            continue;
        }
        let rights: Vec<usize> = (0..parents.len()).filter(|&k| suffix[k] >= n - s).collect();
        if rights.is_empty() {
            continue;
        }
        // A join needs two different parents, so one shared index is not enough.
        if lefts.len() > 1 || rights.len() > 1 || lefts[0] != rights[0] {
            return true;
        }
    }
    false
}

/// Flag every bimera in a dereplicated set, in the order given.
///
/// A sequence is tested against those before it that clear both abundance
/// rules, so the walk has to be in decreasing abundance. The tests are
/// independent once the ordering is fixed, which is what makes this parallel.
#[must_use]
pub fn flag_bimeras(uniques: &[Unique], options: &BimeraOptions) -> Vec<bool> {
    (0..uniques.len())
        .into_par_iter()
        .map(|i| {
            let query = &uniques[i];
            let threshold =
                (query.abundance as f64 * options.min_fold_parent_over_abundance).ceil() as u32;
            let parents: Vec<&Unique> = uniques[..i]
                .iter()
                .filter(|p| {
                    p.abundance >= options.min_parent_abundance && p.abundance >= threshold
                })
                .collect();
            is_bimera(&query.sequence, &parents, options.max_shift)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn u(seq: &str, abundance: u32) -> Unique {
        Unique { sequence: seq.as_bytes().to_vec(), abundance }
    }

    #[test]
    fn a_left_right_join_of_two_parents_is_a_bimera() {
        let left = u("AAAACCCC", 100);
        let right = u("GGGGTTTT", 100);
        let chimera = "AAAATTTT";
        assert!(is_bimera(chimera.as_bytes(), &[&left, &right], 0));
    }

    #[test]
    fn a_sequence_unlike_its_parents_is_not_a_bimera() {
        let left = u("AAAACCCC", 100);
        let right = u("GGGGTTTT", 100);
        assert!(!is_bimera(b"ACACACAC", &[&left, &right], 0));
    }

    #[test]
    fn one_parent_cannot_supply_both_sides() {
        // A sequence identical to a single parent is that parent, not a
        // chimera of it with itself.
        let only = u("AAAACCCC", 100);
        let other = u("AAAACCCC", 90);
        assert!(!is_bimera(b"AAAACCCC", &[&only, &other], 0));
    }

    #[test]
    fn parents_must_clear_the_abundance_rules() {
        // The join is there, but the parents are too rare to be believed.
        let uniques = vec![
            u("AAAACCCC", 4),
            u("GGGGTTTT", 4),
            u("AAAATTTT", 2),
        ];
        let flags = flag_bimeras(&uniques, &BimeraOptions::default());
        assert!(!flags[2], "parents below min_parent_abundance should not convict");

        let abundant = vec![
            u("AAAACCCC", 40),
            u("GGGGTTTT", 40),
            u("AAAATTTT", 2),
        ];
        let flags = flag_bimeras(&abundant, &BimeraOptions::default());
        assert!(flags[2], "abundant parents should convict");
    }

    #[test]
    fn a_shift_lets_parents_of_a_different_length_still_convict() {
        // The left parent has one extra base at its start, so no
        // zero-shift comparison finds the join and a shifted one does.
        let left = u("TAAAACCCC", 100);
        let right = u("GGGGTTTT", 100);
        let chimera = b"AAAATTTT";
        assert!(!is_bimera(chimera, &[&left, &right], 0));
        assert!(is_bimera(chimera, &[&left, &right], 16));
    }

    #[test]
    fn the_most_abundant_sequence_is_never_a_bimera() {
        let uniques = vec![u("AAAATTTT", 100), u("AAAACCCC", 50), u("GGGGTTTT", 50)];
        let flags = flag_bimeras(&uniques, &BimeraOptions::default());
        assert!(!flags[0], "nothing precedes the first sequence");
    }
}
