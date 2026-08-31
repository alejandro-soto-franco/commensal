//! Binding against activation for flagellin at TLR5.
//!
//! A flagellin that binds TLR5 tightly need not switch it on. `RhFlaB` from
//! *Roseburia hominis* binds TLR5 more tightly than *Salmonella* `StFliC` and
//! is silent, so affinity ranks the two the wrong way round.
//!
//! # Assembly architecture
//!
//! Signalling wants the 2:2 assembly of the crystal structure (Yoon et al.
//! 2012), and the literature settles where the discrimination sits. Forstneric
//! et al. 2017 find that deleting the D0 domain leaves binding to the TLR5
//! monomer nearly intact while abolishing activation, and assign D0 the job of
//! "crosslinking the two TLR5 receptor monomers and hence stabilizing the
//! functional signaling complex". D0 does not decide whether the complex forms.
//! It decides how long it lasts once it has.
//!
//! That gives three states and two races:
//!
//! 1. a flagellin binds one receptor, giving a 1:1 complex, at `k_on` and
//!    `k_off`, which is what surface plasmon resonance against the truncated
//!    ectodomain measures;
//! 2. two 1:1 complexes meet in the membrane and D0 crosslinks them into the
//!    2:2 assembly, which is a race against the 1:1 complex falling apart;
//! 3. the assembly recruits its adaptors, which is a race against the 2:2
//!    complex falling apart, and this is the race D0 governs.
//!
//! Because the species that meet in step 2 are both *occupied* receptors, the
//! assembly rate goes as the square of occupancy rather than linearly in
//! receptor density. That is a prediction about the shape of a dose response,
//! not a detail of bookkeeping.
//!
//! # Purpose of the model
//!
//! Every published measurement of these flagellins covers step 1. Nothing
//! measures step 3, and step 3 is where silence is decided. The model takes the
//! measured step-1 kinetics and a measured activation ratio and returns what
//! step 3 must be doing to account for the difference, which is the number a
//! binding assay cannot report.
//!
//! Brownian dynamics enters through the survival distributions. A ligand
//! leaving a pocket does not leave at a single rate: partial escape and
//! rebinding put weight at long times, and a race is decided by the whole
//! distribution rather than its mean. `ermak` supplies that shape and the
//! measured rate sets its timescale.

mod activation;
mod encounter;

pub use activation::{
    Construct, Decomposition, activation_probability, activation_probability_sampled, decompose,
    measured_constructs, required_dimer_lifetime, signalling_probability, survival_sample,
};
pub use encounter::{KT_310K, assembly_rate, encounter_rate, occupancy};
