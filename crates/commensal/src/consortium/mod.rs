//! Colonisation and butyrate production for a defined gut consortium.
//!
//! A synbiotic is a bacterium and the substrate it needs, given together
//! because neither works alone: the strain without its substrate cannot
//! withstand washout, and the substrate without the strain is fermented by
//! whatever is already there. Hesser et al. 2024 give *Anaerostipes caccae*
//! with lactulose and raise luminal butyrate in gnotobiotic mice colonised from
//! an allergic infant, where either component alone does much less.
//!
//! A curve fitted to butyrate against time describes that. It does not say why
//! the combination beats its parts, what happens at a dose nobody ran, or
//! whether the strain persists once dosing stops. Those are questions about a
//! mechanism, so this is the mechanism.
//!
//! # Chemostat model
//!
//! The colon as a chemostat. Substrate arrives, the strain grows on it by Monod
//! kinetics, butyrate comes off the fermentation, and everything is diluted by
//! transit while butyrate is additionally absorbed across the epithelium.
//!
//! ```text
//! dS/dt = I(t) - (mu_max / yield) B S / (Ks + S) - D S
//! dB/dt =        mu_max         B S / (Ks + S) - D B
//! dP/dt = alpha  (mu_max / yield) B S / (Ks + S) - (k_abs + D) P
//! ```
//!
//! `S` is substrate, `B` is the strain, `P` is luminal butyrate. Everything is
//! per unit luminal volume and time is in hours.
//!
//! # Two thresholds the model states
//!
//! Growth cannot outrun transit unless `mu_max > D`, and even then the strain
//! needs the substrate kept above `Ks D / (mu_max - D)` to persist. Below
//! either, the strain washes out however much of it is given, which is the
//! quantitative form of why a probiotic alone fails.

mod chemostat;
mod monod;

pub use chemostat::{Arm, State, simulate, synergy};
pub use monod::{Consortium, SteadyState};
