//! Reads to unique sequences.
//!
//! Two steps of the 16S pipeline, both of which dominate its runtime at scale
//! and neither of which needs any of the statistics around them.
//!
//! Dereplication collapses reads to unique sequences with counts. It is hashing
//! and nothing else, and the reference implementation does it inside R.
//!
//! Bimera detection removes the chimeric sequences that PCR makes when a
//! partial extension from one template primes on another. A sequence is a
//! bimera when it can be read as the left of one more abundant sequence
//! followed by the right of another, so the test is over pairs of parents and
//! over split points, which is quadratic in the number of unique sequences and
//! is where the time goes.

mod bimera;
mod derep;
mod fastq;

pub use bimera::{BimeraOptions, flag_bimeras, is_bimera};
pub use derep::dereplicate;
pub use fastq::read_fastq;
