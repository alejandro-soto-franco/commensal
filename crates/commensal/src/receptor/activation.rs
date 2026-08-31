use ermak::kinetics::escape_times;

/// One flagellin construct as measured by surface plasmon resonance against
/// TLR5N14, the truncated ectodomain presenting the D1 sites alone.
#[derive(Clone, Debug)]
pub struct Construct {
    pub name: &'static str,
    /// Association rate, per molar per second.
    pub k_on: f64,
    /// Dissociation rate, per second.
    pub k_off: f64,
    /// Whether the construct activates TLR5 in a reporter assay.
    pub stimulatory: bool,
}

impl Construct {
    /// The equilibrium dissociation constant in molar, `k_off / k_on`.
    #[must_use]
    pub fn k_d(&self) -> f64 {
        self.k_off / self.k_on
    }

    /// Mean survival of the first complex, in seconds.
    #[must_use]
    pub fn residence_s(&self) -> f64 {
        1.0 / self.k_off
    }
}

/// The four constructs of Bell et al. 2024, measured against TLR5N14.
///
/// The first two are the point of the whole exercise: the silent flagellin has
/// the smaller dissociation constant and the shorter-lived complex.
#[must_use]
pub fn measured_constructs() -> Vec<Construct> {
    vec![
        Construct { name: "StFliC", k_on: 11_293.0, k_off: 1.925e-4, stimulatory: true },
        Construct { name: "RhFlaB", k_on: 45_157.0, k_off: 5.547e-4, stimulatory: false },
        Construct { name: "StFliC-PIM", k_on: 18_830.0, k_off: 1.100e-2, stimulatory: false },
        Construct { name: "RhFlaB-PIM", k_on: 9_091.0, k_off: 6.907e-4, stimulatory: false },
    ]
}

/// Activation probability per binding event, for exponential survival.
///
/// The closed form: the second receptor is engaged at rate `p * k_enc` and the
/// complex dies at `k_off`, so the first of the two decides the outcome.
#[must_use]
pub fn activation_probability(engagement: f64, k_enc: f64, k_off: f64) -> f64 {
    let productive = engagement * k_enc;
    productive / (productive + k_off)
}

/// Activation probability per binding event, over a supplied survival sample.
///
/// `survival` is a sample of complex lifetimes in seconds. Each one gives
/// `1 - exp(-p k_enc T)`, the chance that at least one productive engagement
/// lands inside it, and the average over the sample is the activation
/// probability. This is the form to use when survival is not exponential.
#[must_use]
pub fn activation_probability_sampled(engagement: f64, k_enc: f64, survival: &[f64]) -> f64 {
    assert!(!survival.is_empty(), "an empty survival sample decides nothing");
    let productive = engagement * k_enc;
    survival
        .iter()
        .map(|t| 1.0 - (-productive * t).exp())
        .sum::<f64>()
        / survival.len() as f64
}

/// A Brownian survival sample rescaled to a measured dissociation rate.
///
/// Brownian dynamics in a pocket gives the *shape* of the survival
/// distribution, which is not exponential, and the measurement gives its
/// timescale. Rescaling the sample so its mean is `1 / k_off` keeps both: the
/// simulated non-exponentiality and the measured rate.
#[must_use]
pub fn survival_sample(k_off: f64, barrier: f64, replicas: usize, seed: u64) -> Vec<f64> {
    let raw = escape_times(barrier, 1.0, 1.0, 1e-3, 400_000, replicas, seed);
    let mean = raw.iter().sum::<f64>() / raw.len() as f64;
    let scale = (1.0 / k_off) / mean;
    raw.into_iter().map(|t| t * scale).collect()
}

/// How activation splits between the two lifetimes, for one pair of constructs.
#[derive(Clone, Debug)]
pub struct Decomposition {
    /// The measured ratio of activation between the two.
    pub observed_fold: f64,
    /// What the measured difference in 1:1 lifetime supplies, if activation
    /// depends on it through a single race.
    pub monomer_fold: f64,
    /// What the 2:2 complex must then supply, which is the quantity D0 governs
    /// and which no binding assay reports.
    pub dimer_fold: f64,
    /// The number of sequential races that would let the measured 1:1 kinetics
    /// account for the whole ratio on their own.
    pub proofreading_steps: f64,
}

/// Split an observed activation ratio into a 1:1 part and a 2:2 part.
///
/// Where activation is rare per binding event each race contributes its own
/// lifetime ratio and the two multiply, so the measured `k_off` values fix the
/// first factor and leave the second as the only unknown.
///
/// `proofreading_steps` answers the obvious objection. If activation ran
/// through `n` sequential races against the *same* 1:1 dissociation, the
/// measured ratio would be raised to the `n`th power, and a large enough `n`
/// would explain the observation with no help from D0 at all. The field reports
/// this exponent for the T cell receptor (McKeithan 1995). Reporting the `n`
/// that would be required makes the alternative explicit and testable rather
/// than leaving it implied.
#[must_use]
pub fn decompose(observed_fold: f64, stimulatory: &Construct, silent: &Construct) -> Decomposition {
    let monomer_fold = silent.k_off / stimulatory.k_off;
    Decomposition {
        observed_fold,
        monomer_fold,
        dimer_fold: observed_fold / monomer_fold,
        proofreading_steps: observed_fold.ln() / monomer_fold.ln(),
    }
}

/// The 2:2 complex lifetime a construct needs to reach a given signalling
/// probability, in seconds, when adaptors are recruited at `k_sig`.
///
/// Inverts `P = k_sig / (k_sig + 1 / tau)`, the race of step 3.
#[must_use]
pub fn required_dimer_lifetime(signalling_probability: f64, k_sig: f64) -> f64 {
    assert!(
        (0.0..1.0).contains(&signalling_probability),
        "a probability below one is required"
    );
    signalling_probability / (k_sig * (1.0 - signalling_probability))
}

/// The probability that an assembled 2:2 complex signals before it falls apart.
#[must_use]
pub fn signalling_probability(dimer_lifetime: f64, k_sig: f64) -> f64 {
    let decay = 1.0 / dimer_lifetime;
    k_sig / (k_sig + decay)
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::encounter::encounter_rate;

    #[test]
    fn affinity_ranks_the_silent_flagellin_above_the_agonist() {
        // The fact the whole model exists to explain. Ranking by how tightly a
        // flagellin binds puts the silent one first, so binding alone cannot be
        // what decides activation.
        let c = measured_constructs();
        let (stimulatory, silent) = (&c[0], &c[1]);
        assert!(stimulatory.stimulatory && !silent.stimulatory);
        assert!(
            silent.k_d() < stimulatory.k_d(),
            "the silent flagellin should bind more tightly: {:.2} nM against {:.2} nM",
            silent.k_d() * 1e9,
            stimulatory.k_d() * 1e9
        );
        assert!(
            silent.residence_s() < stimulatory.residence_s(),
            "and should still form the shorter-lived complex"
        );
    }

    #[test]
    fn the_sampled_form_reproduces_the_closed_form_on_exponential_survival() {
        // Draw exponential lifetimes and confirm the average of
        // 1 - exp(-p k_enc T) lands on p k_enc / (p k_enc + k_off).
        let (engagement, k_enc, k_off) = (1e-3, 3.5, 5e-4);
        let mut state = 0x2545_F491_4F6C_DD1D_u64;
        let survival: Vec<f64> = (0..200_000)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                let u = ((state >> 11) as f64 / (1_u64 << 53) as f64).max(1e-16);
                -u.ln() / k_off
            })
            .collect();
        let sampled = activation_probability_sampled(engagement, k_enc, &survival);
        let closed = activation_probability(engagement, k_enc, k_off);
        assert!(
            (sampled - closed).abs() < 5e-3 * closed,
            "sampled {sampled:.6} against closed {closed:.6}"
        );
    }

    #[test]
    fn a_longer_lived_complex_activates_more_often() {
        let k_enc = encounter_rate(0.2, 10.0, 0.005);
        let slow = activation_probability(1e-4, k_enc, 1.925e-4);
        let fast = activation_probability(1e-4, k_enc, 5.547e-4);
        assert!(slow > fast, "{slow} should exceed {fast}");
    }

    #[test]
    fn the_decomposition_multiplies_back_to_what_was_observed() {
        let c = measured_constructs();
        let d = decompose(100.0, &c[0], &c[1]);
        assert!((d.monomer_fold * d.dimer_fold - 100.0).abs() < 1e-9);
        assert!(
            d.monomer_fold > 2.5 && d.monomer_fold < 3.5,
            "the 1:1 lifetime ratio is about threefold, found {}",
            d.monomer_fold
        );
        // And the alternative, stated so it can be argued with.
        assert!(
            d.proofreading_steps > 4.0 && d.proofreading_steps < 5.0,
            "about four sequential races would be needed instead, found {}",
            d.proofreading_steps
        );
    }

    #[test]
    fn the_silent_construct_binds_tighter_and_unbinds_sooner() {
        // The result the module exists to state. The silent flagellin has the
        // higher affinity, 12.3 nM against 17.1 nM, and it reaches that
        // through an on-rate four times faster while unbinding nearly three
        // times sooner. Weak binding therefore cannot explain the silence,
        // which is what sends the model to complex lifetime instead.
        let constructs = measured_constructs();
        let silent = constructs.iter().find(|c| !c.stimulatory).expect("a silent construct");
        let stim = constructs.iter().find(|c| c.stimulatory).expect("a stimulatory construct");
        assert!((silent.k_d() * 1e9 - 12.3).abs() < 0.1, "silent K_d {}", silent.k_d() * 1e9);
        assert!((stim.k_d() * 1e9 - 17.1).abs() < 0.1, "stimulatory K_d {}", stim.k_d() * 1e9);
        let ratio = silent.residence_s() / stim.residence_s();
        assert!((ratio - 0.347).abs() < 0.01, "lifetime ratio {ratio}");
        assert!(
            silent.residence_s() < stim.residence_s(),
            "the silent construct should unbind sooner: {} against {}",
            silent.residence_s(),
            stim.residence_s()
        );
        assert!(
            (silent.k_on / stim.k_on - 4.0).abs() < 0.01,
            "on-rate ratio {}",
            silent.k_on / stim.k_on
        );
    }
}
