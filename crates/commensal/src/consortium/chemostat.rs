use super::monod::Consortium;

/// The luminal state at one instant.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct State {
    pub substrate: f64,
    pub biomass: f64,
    pub butyrate: f64,
}

/// One trajectory under a substrate supply that may vary with time.
///
/// Integrated with classical Runge-Kutta at fixed step. `supply` is the
/// substrate arriving per unit volume per hour, which is how a gavage schedule
/// or a continuous dietary source enters.
#[must_use]
pub fn simulate<F>(
    consortium: &Consortium,
    initial: State,
    supply: F,
    hours: f64,
    dt: f64,
) -> Vec<(f64, State)>
where
    F: Fn(f64) -> f64,
{
    assert!(dt > 0.0 && hours > 0.0);
    let derivative = |t: f64, s: State| -> State {
        let uptake = consortium.growth_rate(s.substrate) * s.biomass
            / consortium.yield_coefficient;
        State {
            substrate: supply(t) - uptake - consortium.dilution * s.substrate,
            biomass: consortium.growth_rate(s.substrate) * s.biomass
                - consortium.dilution * s.biomass,
            butyrate: consortium.butyrate_per_substrate * uptake
                - (consortium.absorption + consortium.dilution) * s.butyrate,
        }
    };
    let step = |t: f64, s: State, h: f64| -> State {
        let add = |a: State, b: State, w: f64| State {
            substrate: a.substrate + w * b.substrate,
            biomass: a.biomass + w * b.biomass,
            butyrate: a.butyrate + w * b.butyrate,
        };
        let k1 = derivative(t, s);
        let k2 = derivative(t + h / 2.0, add(s, k1, h / 2.0));
        let k3 = derivative(t + h / 2.0, add(s, k2, h / 2.0));
        let k4 = derivative(t + h, add(s, k3, h));
        State {
            substrate: (s.substrate
                + h / 6.0 * (k1.substrate + 2.0 * k2.substrate + 2.0 * k3.substrate + k4.substrate))
                .max(0.0),
            biomass: (s.biomass
                + h / 6.0 * (k1.biomass + 2.0 * k2.biomass + 2.0 * k3.biomass + k4.biomass))
                .max(0.0),
            butyrate: (s.butyrate
                + h / 6.0 * (k1.butyrate + 2.0 * k2.butyrate + 2.0 * k3.butyrate + k4.butyrate))
                .max(0.0),
        }
    };

    let steps = (hours / dt).ceil() as usize;
    let mut out = Vec::with_capacity(steps + 1);
    let mut s = initial;
    out.push((0.0, s));
    for i in 0..steps {
        let t = i as f64 * dt;
        s = step(t, s, dt);
        out.push((t + dt, s));
    }
    out
}

/// The four arms of a synbiotic experiment.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Arm {
    /// Neither component.
    Untreated,
    /// The strain alone.
    Probiotic,
    /// The substrate alone.
    Prebiotic,
    /// Both.
    Synbiotic,
}

/// Peak luminal butyrate reached in each arm, and the excess of the combination
/// over what its parts predict.
///
/// The excess is what a curve fit cannot produce. Two components that each do
/// nothing can together do a great deal, and the size of that gap is a property
/// of the mechanism rather than of the data.
#[must_use]
pub fn synergy(
    consortium: &Consortium,
    dose: f64,
    inoculum: f64,
    resident_butyrate: f64,
    hours: f64,
) -> (Vec<(Arm, f64)>, f64) {
    let peak = |arm: Arm| -> f64 {
        let biomass = match arm {
            Arm::Probiotic | Arm::Synbiotic => inoculum,
            _ => 0.0,
        };
        let supply = move |_: f64| match arm {
            Arm::Prebiotic | Arm::Synbiotic => dose,
            _ => 0.0,
        };
        let initial = State { substrate: 0.0, biomass, butyrate: resident_butyrate };
        simulate(consortium, initial, supply, hours, 0.01)
            .into_iter()
            .map(|(_, s)| s.butyrate)
            .fold(0.0_f64, f64::max)
    };

    let arms: Vec<(Arm, f64)> = [Arm::Untreated, Arm::Probiotic, Arm::Prebiotic, Arm::Synbiotic]
        .into_iter()
        .map(|a| (a, peak(a)))
        .collect();

    let base = arms[0].1;
    let alone = (arms[1].1 - base) + (arms[2].1 - base);
    let together = arms[3].1 - base;
    let excess = if alone.abs() < 1e-12 {
        f64::INFINITY
    } else {
        together / alone
    };
    (arms, excess)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * b.abs().max(1.0)
    }

    #[test]
    fn the_simulation_settles_on_the_steady_state_it_predicts() {
        // The closed form and the integrator have to agree, or one of them is
        // wrong. This is the test that pins both.
        let c = Consortium::default();
        let input = 0.4;
        let predicted = c.steady_state(input);
        assert!(predicted.persists);

        let trace = simulate(
            &c,
            State { substrate: 0.1, biomass: 0.02, butyrate: 0.0 },
            |_| input,
            4000.0,
            0.01,
        );
        let (_, last) = trace.last().copied().expect("a trajectory has an end");
        assert!(
            approx(last.substrate, predicted.substrate, 1e-3),
            "substrate {} against {}",
            last.substrate,
            predicted.substrate
        );
        assert!(
            approx(last.biomass, predicted.biomass, 1e-3),
            "biomass {} against {}",
            last.biomass,
            predicted.biomass
        );
        assert!(
            approx(last.butyrate, predicted.butyrate, 1e-3),
            "butyrate {} against {}",
            last.butyrate,
            predicted.butyrate
        );
    }

    #[test]
    fn the_strain_washes_out_once_dosing_stops() {
        let c = Consortium::default();
        let trace = simulate(
            &c,
            State { substrate: 0.0, biomass: 0.5, butyrate: 0.0 },
            |t| if t < 100.0 { 0.5 } else { 0.0 },
            600.0,
            0.01,
        );
        let during = trace.iter().find(|(t, _)| *t >= 99.0).expect("during").1;
        let after = trace.last().expect("after").1;
        assert!(
            after.biomass < 0.05 * during.biomass,
            "biomass should collapse: {} against {}",
            after.biomass,
            during.biomass
        );
    }

    #[test]
    fn the_combination_beats_the_sum_of_its_parts() {
        let c = Consortium::default();
        let (arms, excess) = synergy(&c, 0.5, 0.02, 0.0, 400.0);
        let peak = |a: Arm| arms.iter().find(|(x, _)| *x == a).expect("arm").1;
        assert!(peak(Arm::Probiotic) < 1e-6, "the strain alone has nothing to ferment");
        assert!(peak(Arm::Prebiotic) < 1e-6, "the substrate alone has nobody to ferment it");
        assert!(peak(Arm::Synbiotic) > 0.05, "together they produce butyrate");
        assert!(
            excess > 100.0 || excess.is_infinite(),
            "the excess over the parts should be large, found {excess}"
        );
    }
}
