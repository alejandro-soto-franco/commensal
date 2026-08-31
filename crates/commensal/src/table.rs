//! The objects more than one module needs.

use std::fmt;

/// Anything that can go wrong in this crate.
#[derive(Debug)]
pub enum Error {
    /// The count buffer does not match the sample and sequence counts.
    Shape { samples: usize, sequences: usize, counts: usize },
    /// Every sample was empty, so there is nothing to analyse.
    Empty,
    /// A file could not be read or written.
    Io(std::io::Error),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Shape { samples, sequences, counts } => write!(
                f,
                "{samples} samples by {sequences} sequences needs {} counts, given {counts}",
                samples * sequences
            ),
            Error::Empty => write!(f, "every sample is empty"),
            Error::Io(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for Error {}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Error::Io(e)
    }
}

/// A unique sequence and how many reads collapsed onto it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Unique {
    pub sequence: Vec<u8>,
    pub abundance: u32,
}

/// A sample by sequence count table.
///
/// Construction refuses any sample whose row sums to zero. Such a sample has no
/// composition: a relative abundance divides by its total, which gives NaN, and
/// a permutation test over a distance matrix containing NaN reports a small
/// p-value that means nothing. Refusing it here defends every consumer at once.
#[derive(Clone, Debug)]
pub struct SequenceTable {
    samples: Vec<String>,
    sequences: Vec<Vec<u8>>,
    counts: Vec<u32>,
}

impl SequenceTable {
    /// Build a table, dropping empty samples and reporting how many went.
    ///
    /// `counts` is row major, one row per sample.
    ///
    /// # Errors
    /// [`Error::Shape`] when `counts` does not have one entry per sample and
    /// sequence, and [`Error::Empty`] when no sample survives.
    pub fn new(
        samples: Vec<String>,
        sequences: Vec<Vec<u8>>,
        counts: Vec<u32>,
    ) -> Result<(Self, usize), Error> {
        let (n, m) = (samples.len(), sequences.len());
        if counts.len() != n * m {
            return Err(Error::Shape { samples: n, sequences: m, counts: counts.len() });
        }
        let live: Vec<usize> = (0..n)
            .filter(|i| counts[i * m..(i + 1) * m].iter().any(|c| *c > 0))
            .collect();
        if live.is_empty() {
            return Err(Error::Empty);
        }
        let dropped = n - live.len();
        let kept_samples = live.iter().map(|i| samples[*i].clone()).collect();
        let mut kept = Vec::with_capacity(live.len() * m);
        for i in &live {
            kept.extend_from_slice(&counts[i * m..(i + 1) * m]);
        }
        Ok((Self { samples: kept_samples, sequences, counts: kept }, dropped))
    }

    #[must_use]
    pub fn samples(&self) -> &[String] {
        &self.samples
    }

    #[must_use]
    pub fn sequences(&self) -> &[Vec<u8>] {
        &self.sequences
    }

    /// The counts, row major.
    #[must_use]
    pub fn counts(&self) -> &[u32] {
        &self.counts
    }

    /// Samples by sequences.
    #[must_use]
    pub fn shape(&self) -> (usize, usize) {
        (self.samples.len(), self.sequences.len())
    }

    /// One sample's counts.
    #[must_use]
    pub fn row(&self, i: usize) -> &[u32] {
        let m = self.sequences.len();
        &self.counts[i * m..(i + 1) * m]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_sample_with_no_reads_is_refused_and_counted() {
        // Filtering can empty a sample. Such a sample has no composition, and
        // a distance that divides by its total returns NaN, which a
        // permutation test then reads as a significant result out of nothing.
        let samples = vec!["a".to_string(), "empty".to_string(), "b".to_string()];
        let sequences = vec![b"ACGT".to_vec(), b"TTTT".to_vec()];
        let counts = vec![3, 1, 0, 0, 2, 5];
        let (table, dropped) =
            SequenceTable::new(samples, sequences, counts).expect("two live samples");
        assert_eq!(dropped, 1);
        assert_eq!(table.shape(), (2, 2));
        assert_eq!(table.samples(), ["a", "b"]);
        assert_eq!(table.row(0), [3, 1]);
        assert_eq!(table.row(1), [2, 5]);
    }

    #[test]
    fn a_table_with_no_live_sample_is_an_error() {
        let samples = vec!["a".to_string()];
        let sequences = vec![b"ACGT".to_vec()];
        assert!(SequenceTable::new(samples, sequences, vec![0]).is_err());
    }

    #[test]
    fn a_count_length_that_does_not_match_the_shape_is_an_error() {
        let samples = vec!["a".to_string()];
        let sequences = vec![b"ACGT".to_vec(), b"TTTT".to_vec()];
        assert!(SequenceTable::new(samples, sequences, vec![1]).is_err());
    }
}
