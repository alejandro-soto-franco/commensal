/// Boltzmann's constant times body temperature, in kcal/mol.
pub const KT_310K: f64 = 0.6161;

/// The rate at which a bound complex meets a second receptor in the membrane,
/// per second.
///
/// Two-dimensional diffusion has no steady-state Smoluchowski rate, so the
/// standard finite-outer-boundary form is used with the boundary at the mean
/// half-spacing of the receptors, `b = 1 / sqrt(pi * density)`:
///
/// ```text
/// k_enc = 2 pi D density / ln(b / a)
/// ```
///
/// `d_lateral` is the relative diffusion of the two receptors in the membrane
/// (twice a single receptor's), in square micrometres per second, `density` is
/// receptors per square micrometre, and `encounter_radius` is in micrometres.
#[must_use]
pub fn encounter_rate(d_lateral: f64, density: f64, encounter_radius: f64) -> f64 {
    assert!(density > 0.0 && encounter_radius > 0.0 && d_lateral > 0.0);
    let outer = 1.0 / (std::f64::consts::PI * density).sqrt();
    let ratio = (outer / encounter_radius).max(std::f64::consts::E.sqrt());
    2.0 * std::f64::consts::PI * d_lateral * density / ratio.ln()
}

/// The fraction of receptors sitting in a 1:1 complex at equilibrium.
///
/// `ligand` and `k_d` share units. Only occupied receptors can assemble, so
/// this is what the assembly step counts rather than the total density.
#[must_use]
pub fn occupancy(ligand: f64, k_d: f64) -> f64 {
    ligand / (ligand + k_d)
}

/// The rate at which one 1:1 complex meets another, per second.
///
/// The physiological assembly is 2:2, so both partners in the encounter are
/// occupied receptors. The rate therefore goes as the square of occupancy, and
/// a dose response that is quadratic at low ligand is a signature of this
/// architecture rather than of cooperativity in binding.
#[must_use]
pub fn assembly_rate(
    d_lateral: f64,
    receptor_density: f64,
    occupancy: f64,
    encounter_radius: f64,
) -> f64 {
    let occupied = receptor_density * occupancy;
    if occupied <= 0.0 {
        return 0.0;
    }
    encounter_rate(d_lateral, occupied, encounter_radius) * occupancy
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_encounter_rate_grows_with_density_and_diffusion() {
        let base = encounter_rate(0.2, 10.0, 0.005);
        assert!(encounter_rate(0.2, 40.0, 0.005) > base);
        assert!(encounter_rate(0.8, 10.0, 0.005) > base);
    }
}
