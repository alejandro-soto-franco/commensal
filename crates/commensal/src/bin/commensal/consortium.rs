//! What the synbiotic does, and what a curve fitted to it cannot say.

use commensal::consortium::{Arm, Consortium, State, simulate, synergy};

pub fn run() {
    let c = Consortium::default();

    println!("A. caccae on lactulose in the colon\n");
    println!("  mu_max            {:>8.3} /h", c.mu_max);
    println!("  dilution          {:>8.3} /h", c.dilution);
    println!("  absorption        {:>8.3} /h", c.absorption);
    match c.washout_substrate() {
        Some(s) => println!(
            "\n  Growth outruns transit, but only above {s:.4} of substrate. Below that the\n  \
             strain washes out however much is given, which is why a probiotic on its own\n  \
             does not hold."
        ),
        None => println!("\n  Growth cannot outrun transit at any substrate concentration."),
    }

    // The four arms.
    let (arms, excess) = synergy(&c, 0.5, 0.02, 0.0, 400.0);
    println!("\nPeak luminal butyrate by arm");
    for (arm, peak) in &arms {
        let name = match arm {
            Arm::Untreated => "untreated",
            Arm::Probiotic => "strain alone",
            Arm::Prebiotic => "substrate alone",
            Arm::Synbiotic => "synbiotic",
        };
        println!("  {name:<18} {peak:>10.4}");
    }
    if excess.is_infinite() {
        println!(
            "\n  Neither component alone raises butyrate at all, so the combination is not a\n  \
             sum of effects. There is nothing for an additive model to add."
        );
    } else {
        println!("\n  The combination exceeds the sum of its parts {excess:.0}-fold.");
    }

    // Dose response, and where it stops paying.
    println!("\nSteady state against lactulose supply");
    println!(
        "{:>10} {:>12} {:>12} {:>12} {:>10}",
        "supply", "substrate", "biomass", "butyrate", "persists"
    );
    for input in [0.005_f64, 0.01, 0.05, 0.2, 0.5, 1.0, 2.0] {
        let s = c.steady_state(input);
        println!(
            "{:>10.3} {:>12.4} {:>12.4} {:>12.4} {:>10}",
            input,
            s.substrate,
            s.biomass,
            s.butyrate,
            if s.persists { "yes" } else { "no" }
        );
    }
    println!(
        "\n  Above the threshold the residual substrate does not move: every extra unit of\n  \
         supply goes into biomass and butyrate rather than accumulating. A dose response\n  \
         measured in butyrate is therefore linear in supply, not saturating, until\n  \
         something else becomes limiting. That is a testable claim about her dosing curve."
    );

    // Persistence after dosing stops.
    let trace = simulate(
        &c,
        State { substrate: 0.0, biomass: 0.5, butyrate: 0.0 },
        |t| if t < 168.0 { 0.5 } else { 0.0 },
        500.0,
        0.01,
    );
    let at = |target: f64| {
        trace
            .iter()
            .find(|(t, _)| *t >= target)
            .map(|(_, s)| *s)
            .expect("within the trace")
    };
    println!("\nAfter a week of dosing, then nothing");
    println!("{:>10} {:>12} {:>12}", "hours", "biomass", "butyrate");
    for t in [24.0_f64, 168.0, 192.0, 240.0, 336.0, 480.0] {
        let s = at(t);
        println!("{:>10.0} {:>12.5} {:>12.5}", t, s.biomass, s.butyrate);
    }
    let during = at(167.0);
    let after = at(240.0);
    println!(
        "\n  Butyrate is gone within a day of the last dose, because absorption clears it\n  \
         {:.0} times faster than transit. The strain falls more slowly but it does fall:\n  \
         {:.0}% of its standing biomass is left three days later. Colonisation here is a\n  \
         supply, not an establishment, and stopping the substrate ends it.",
        c.absorption / c.dilution,
        100.0 * after.biomass / during.biomass
    );
}
