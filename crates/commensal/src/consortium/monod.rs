/// Parameters of one strain fermenting one substrate in the lumen.
#[derive(Clone, Debug)]
pub struct Consortium {
    /// Maximum specific growth rate, per hour.
    pub mu_max: f64,
    /// Half-saturation constant for growth, in the units of `S`.
    pub ks: f64,
    /// Biomass produced per unit substrate consumed.
    pub yield_coefficient: f64,
    /// Butyrate produced per unit substrate consumed.
    pub butyrate_per_substrate: f64,
    /// Dilution by transit, per hour.
    pub dilution: f64,
    /// Absorption of butyrate across the epithelium, per hour.
    pub absorption: f64,
}

impl Default for Consortium {
    /// *A. caccae* on lactulose in a mouse colon, at the order of magnitude the
    /// literature supports. Growth rates for gut anaerobes on a preferred sugar
    /// sit in the tenths per hour, mouse colonic transit gives a dilution of
    /// order one per day, and colonocytes clear short chain fatty acids fast
    /// enough that absorption dominates dilution.
    fn default() -> Self {
        Consortium {
            mu_max: 0.35,
            ks: 0.5,
            yield_coefficient: 0.35,
            butyrate_per_substrate: 0.42,
            dilution: 0.05,
            absorption: 1.2,
        }
    }
}

impl Consortium {
    /// The substrate concentration at which growth exactly balances transit.
    ///
    /// Below this the strain washes out no matter how much is given. `None`
    /// when growth cannot outrun transit at any substrate concentration, which
    /// is the harder failure.
    #[must_use]
    pub fn washout_substrate(&self) -> Option<f64> {
        if self.mu_max <= self.dilution {
            return None;
        }
        Some(self.ks * self.dilution / (self.mu_max - self.dilution))
    }

    /// Specific growth rate at a substrate concentration, per hour.
    #[must_use]
    pub fn growth_rate(&self, substrate: f64) -> f64 {
        self.mu_max * substrate / (self.ks + substrate)
    }

    /// The steady state reached under a constant substrate supply `input`.
    ///
    /// Where the strain persists, substrate is pinned at the washout
    /// concentration whatever the supply, and the surplus goes into biomass and
    /// butyrate. That is the chemostat's defining behaviour and it is why the
    /// model predicts a dose response that saturates in substrate but not in
    /// butyrate.
    #[must_use]
    pub fn steady_state(&self, input: f64) -> SteadyState {
        match self.washout_substrate() {
            Some(s_star) if input > self.dilution * s_star => {
                let consumed = input - self.dilution * s_star;
                let biomass = self.yield_coefficient * consumed / self.dilution;
                let butyrate = self.butyrate_per_substrate * consumed
                    / (self.absorption + self.dilution);
                SteadyState { substrate: s_star, biomass, butyrate, persists: true }
            }
            _ => SteadyState {
                substrate: input / self.dilution,
                biomass: 0.0,
                butyrate: 0.0,
                persists: false,
            },
        }
    }
}

/// Where a constant supply leaves the lumen once transients have gone.
#[derive(Clone, Debug, PartialEq)]
pub struct SteadyState {
    pub substrate: f64,
    pub biomass: f64,
    pub butyrate: f64,
    /// Whether the strain holds against washout.
    pub persists: bool,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * b.abs().max(1.0)
    }

    #[test]
    fn a_strain_that_cannot_outgrow_transit_has_no_washout_threshold() {
        let c = Consortium { mu_max: 0.04, dilution: 0.05, ..Default::default() };
        assert!(c.washout_substrate().is_none());
        assert!(!c.steady_state(10.0).persists);
    }

    #[test]
    fn substrate_is_pinned_at_the_washout_concentration_however_much_arrives() {
        // The chemostat result: above the threshold, supply raises biomass and
        // butyrate but leaves the residual substrate where it was.
        let c = Consortium::default();
        let low = c.steady_state(0.3);
        let high = c.steady_state(3.0);
        assert!(low.persists && high.persists);
        assert!(approx(low.substrate, high.substrate, 1e-12));
        assert!(high.biomass > 9.0 * low.biomass);
        assert!(high.butyrate > 9.0 * low.butyrate);
    }
}
