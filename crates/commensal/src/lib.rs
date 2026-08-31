//! Gut commensal analysis: amplicon sequences, consortium dynamics, and host
//! receptor kinetics.
//!
//! The work splits into a sequence half, which says which organisms a sample
//! contains, and a kinetic half, which says what they produce and how the host
//! reads them. They meet at flagellin: a defined consortium has known flagellin
//! sequences, and whether a given one activates TLR5 decides what the host
//! sees.

mod table;

pub mod amplicon;

pub use table::{Error, SequenceTable, Unique};
