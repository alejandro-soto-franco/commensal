//! What binding says, what activation needs, and how far apart the two are.

use commensal::receptor::{
    activation_probability, activation_probability_sampled, assembly_rate, decompose,
    encounter_rate, measured_constructs, occupancy, required_dimer_lifetime, survival_sample,
};

/// The D0 swap of Clasen et al. 2023: giving RhFlaB the D0 domain of FliC makes
/// it a hundredfold more stimulatory, with the D1 binding site untouched.
const D0_SWAP_FOLD: f64 = 100.0;

/// Membrane parameters. TLR5 lateral diffusion is of order 0.1 square
/// micrometres per second, so two receptors approach at twice that, and the
/// encounter radius is a receptor's own size.
const D_LATERAL: f64 = 0.2;
const ENCOUNTER_RADIUS: f64 = 0.005;

pub fn run() {
    let constructs = measured_constructs();

    println!("Measured against TLR5N14, the D1 sites alone\n");
    println!(
        "{:<12} {:>10} {:>12} {:>12} {:>10}",
        "construct", "K_D (nM)", "k_on (1/Ms)", "k_off (1/s)", "1/k_off"
    );
    for c in &constructs {
        println!(
            "{:<12} {:>10.2} {:>12.0} {:>12.3e} {:>8.0} s{}",
            c.name,
            c.k_d() * 1e9,
            c.k_on,
            c.k_off,
            c.residence_s(),
            if c.stimulatory { "  stimulatory" } else { "" }
        );
    }

    let stim = &constructs[0];
    let silent = &constructs[1];
    println!(
        "\nAffinity puts the silent flagellin first: {:.2} nM against {:.2} nM.",
        silent.k_d() * 1e9,
        stim.k_d() * 1e9
    );
    println!(
        "Complex lifetime puts the agonist first: {:.0} s against {:.0} s, a factor {:.2}.",
        stim.residence_s(),
        silent.residence_s(),
        silent.k_off / stim.k_off
    );

    // How often a bound complex meets a second receptor, across the densities a
    // reporter line and an epithelium span.
    println!("\nSecond-receptor encounter rate against receptor density");
    println!("{:>12} {:>14}", "per um^2", "k_enc (1/s)");
    for density in [0.5_f64, 2.0, 10.0, 50.0, 200.0] {
        println!(
            "{:>12.1} {:>14.3}",
            density,
            encounter_rate(D_LATERAL, density, ENCOUNTER_RADIUS)
        );
    }

    let k_enc = encounter_rate(D_LATERAL, 10.0, ENCOUNTER_RADIUS);
    println!(
        "\nAt 10 receptors per um^2 a complex meets one every {:.2} s, so the agonist",
        1.0 / k_enc
    );
    println!(
        "gets about {:.0} encounters before it dissociates and the silent flagellin {:.0}.",
        k_enc / stim.k_off,
        k_enc / silent.k_off
    );

    // The decomposition.
    let d = decompose(D0_SWAP_FOLD, stim, silent);
    println!("\nWhere the hundredfold of the D0 swap comes from");
    println!("  observed                          {:>8.1}x", d.observed_fold);
    println!("  the measured 1:1 lifetime supplies {:>7.2}x", d.monomer_fold);
    println!("  the 2:2 complex must supply       {:>8.1}x", d.dimer_fold);
    println!(
        "\nThe 1:1 kinetics account for {:.0}% of the gap on a log scale. D0 leaves binding\n\
         to the monomer nearly intact, so the rest sits in how long the assembled complex\n\
         lasts, which is the quantity no published measurement covers.",
        100.0 * d.monomer_fold.ln() / d.observed_fold.ln()
    );
    println!(
        "\nThe alternative, stated so it can be tested: {:.2} sequential races against the\n\
         same 1:1 dissociation would account for the whole ratio with no help from D0.\n\
         That is the T cell receptor's proofreading architecture applied to TLR5, and it\n\
         predicts a {:.1}x spread from the PIM series that a single race does not.",
        d.proofreading_steps,
        (constructs[2].k_off / stim.k_off).powf(d.proofreading_steps)
    );

    // What the 2:2 complex has to do, in seconds.
    println!("\nWhat the 2:2 lifetime has to be, for adaptor recruitment at k_sig");
    println!("{:>12} {:>16} {:>16}", "k_sig (1/s)", "StFliC tau (s)", "RhFlaB tau (s)");
    for k_sig in [0.01_f64, 0.1, 1.0] {
        let a = required_dimer_lifetime(0.5, k_sig);
        let b = required_dimer_lifetime(0.5 / d.dimer_fold, k_sig);
        println!("{k_sig:>12.2} {a:>16.2} {b:>16.3}");
    }

    // Brownian survival rather than a single rate.
    println!("\nActivation per binding event, engagement probability 1e-4");
    println!(
        "{:<12} {:>14} {:>16} {:>10}",
        "construct", "closed form", "Brownian sample", "ratio"
    );
    for c in [stim, silent] {
        let survival = survival_sample(c.k_off, 4.0, 4096, 20260831);
        let closed = activation_probability(1e-4, k_enc, c.k_off);
        let sampled = activation_probability_sampled(1e-4, k_enc, &survival);
        println!(
            "{:<12} {:>14.4} {:>16.4} {:>10.3}",
            c.name, closed, sampled, sampled / closed
        );
    }

    // The dose response the 2:2 architecture predicts.
    println!("\nAssembly rate against ligand, which goes as occupancy squared");
    println!("{:>12} {:>12} {:>16} {:>10}", "ligand (nM)", "occupancy", "k_assembly (1/s)", "slope");
    let mut previous: Option<(f64, f64)> = None;
    for ligand_nm in [0.1_f64, 1.0, 10.0, 100.0, 1000.0] {
        let ligand = ligand_nm * 1e-9;
        let theta = occupancy(ligand, stim.k_d());
        let rate = assembly_rate(D_LATERAL, 10.0, theta, ENCOUNTER_RADIUS);
        let slope = previous
            .map(|(l, r): (f64, f64)| (rate / r).ln() / (ligand / l).ln())
            .unwrap_or(f64::NAN);
        println!("{ligand_nm:>12.1} {theta:>12.4} {rate:>16.4} {slope:>10.2}");
        previous = Some((ligand, rate));
    }
    println!(
        "A slope near two below saturation is the signature of the 2:2 architecture, and\n\
         it distinguishes this model from any that assembles on a single receptor."
    );
}
