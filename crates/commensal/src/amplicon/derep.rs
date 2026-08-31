use rustc_hash::FxHashMap;

use crate::Unique;

/// Collapse reads to unique sequences, most abundant first.
///
/// Ties are broken by the sequence itself so the order is total and a run is
/// reproducible. Downstream the order matters: bimera detection walks it, and
/// a sequence may only be called against parents that precede it.
#[must_use]
pub fn dereplicate<I, S>(reads: I) -> Vec<Unique>
where
    I: IntoIterator<Item = S>,
    S: AsRef<[u8]>,
{
    let mut counts: FxHashMap<Vec<u8>, u32> = FxHashMap::default();
    for read in reads {
        *counts.entry(read.as_ref().to_vec()).or_insert(0) += 1;
    }
    let mut uniques: Vec<Unique> = counts
        .into_iter()
        .map(|(sequence, abundance)| Unique { sequence, abundance })
        .collect();
    uniques.sort_unstable_by(|a, b| {
        b.abundance
            .cmp(&a.abundance)
            .then_with(|| a.sequence.cmp(&b.sequence))
    });
    uniques
}

#[cfg(test)]
mod tests {
    use super::*;

    fn u(seq: &str, abundance: u32) -> Unique {
        Unique { sequence: seq.as_bytes().to_vec(), abundance }
    }

    #[test]
    fn dereplication_counts_and_orders() {
        let reads = vec!["ACGT", "ACGT", "TTTT", "ACGT", "GGGG", "TTTT"];
        let d = dereplicate(reads);
        assert_eq!(d[0], u("ACGT", 3));
        assert_eq!(d[1], u("TTTT", 2));
        assert_eq!(d[2], u("GGGG", 1));
    }

    #[test]
    fn ties_break_on_the_sequence_so_a_run_repeats() {
        let a = dereplicate(vec!["TTTT", "AAAA"]);
        let b = dereplicate(vec!["AAAA", "TTTT"]);
        assert_eq!(a, b);
        assert_eq!(a[0].sequence, b"AAAA".to_vec());
    }
}
