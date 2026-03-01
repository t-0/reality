use std::f64::consts::PI;

use crate::couplings::constants::*;
use crate::couplings::running::alpha_s;
use crate::particles::PdgId;

use super::amplitude::MatrixElement;

/// e+e- -> q qbar at tree level.
///
/// This is the same as e+e- -> mu+mu- but with:
/// - Color factor N_c = 3
/// - Quark charge Q_q
/// - Optional O(alpha_s) QCD correction factor: (1 + alpha_s/pi)
pub struct EpEmToQQbar {
    pub quark_pdg: PdgId,
    pub quark_charge: f64,
    pub quark_mass: f64,
    pub include_qcd_correction: bool,
}

impl EpEmToQQbar {
    pub fn new(quark_pdg: PdgId, charge: f64, mass: f64) -> Self {
        Self {
            quark_pdg,
            quark_charge: charge,
            quark_mass: mass,
            include_qcd_correction: true,
        }
    }

    /// Up-type quark (u, c): charge = +2/3
    pub fn up() -> Self {
        Self::new(PdgId::UP, 2.0 / 3.0, 0.00216)
    }

    pub fn down() -> Self {
        Self::new(PdgId::DOWN, -1.0 / 3.0, 0.00467)
    }

    pub fn charm() -> Self {
        Self::new(PdgId::CHARM, 2.0 / 3.0, 1.27)
    }

    pub fn strange() -> Self {
        Self::new(PdgId::STRANGE, -1.0 / 3.0, 0.0934)
    }

    pub fn bottom() -> Self {
        Self::new(PdgId::BOTTOM, -1.0 / 3.0, 4.18)
    }
}

impl MatrixElement for EpEmToQQbar {
    fn amplitude_squared(&self, s: f64, t: f64) -> f64 {
        // Same as e+e- -> mu+mu- but with Q_f = quark_charge and N_c = 3
        let u = -s - t;
        let e4 = (4.0 * PI * ALPHA_EM).powi(2);

        // |M|^2_avg for e+e- -> f fbar (QED) = 2 * e^4 * Q_f^4 * (t^2+u^2)/s^2
        // Wait: the charge factors are Q_e^2 * Q_f^2 in the propagator.
        // Since Q_e = -1, Q_e^2 = 1.
        // |M|^2_avg = 2 * e^4 * Q_q^2 * (t^2 + u^2) / s^2
        let q2 = self.quark_charge * self.quark_charge;
        let mut m2 = N_C * 2.0 * e4 * q2 * (t * t + u * u) / (s * s);

        // O(alpha_s) QCD correction
        if self.include_qcd_correction {
            let as_val = alpha_s(s.sqrt());
            m2 *= 1.0 + as_val / PI;
        }

        m2
    }

    fn initial_state(&self) -> Vec<PdgId> {
        vec![PdgId::POSITRON, PdgId::ELECTRON]
    }

    fn final_state(&self) -> Vec<PdgId> {
        vec![self.quark_pdg.antiparticle(), self.quark_pdg]
    }

    fn initial_state_dof(&self) -> f64 {
        4.0 // 2 spins x 2 spins (initial e+e-)
    }

    fn final_masses(&self) -> (f64, f64) {
        (self.quark_mass, self.quark_mass)
    }

    fn name(&self) -> &str {
        "e+e- -> q qbar"
    }
}

/// Compute the R-ratio at a given sqrt(s):
/// R = sigma(e+e- -> hadrons) / sigma(e+e- -> mu+mu-)
///   = N_c * sum_q Q_q^2 (for quarks with 2*m_q < sqrt_s)
pub fn r_ratio(sqrt_s: f64) -> f64 {
    let charges = [
        (2.0 / 3.0, 0.00216),  // u
        (-1.0 / 3.0, 0.00467), // d
        (-1.0 / 3.0, 0.0934),  // s
        (2.0 / 3.0, 1.27),     // c
        (-1.0 / 3.0, 4.18),    // b
        (2.0 / 3.0, 172.57),   // t
    ];

    let mut r = 0.0;
    for &(q, m) in &charges {
        if sqrt_s > 2.0 * m {
            r += q * q;
        }
    }
    N_C * r
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_r_ratio_below_charm() {
        // Below charm threshold (~2.5 GeV), only u, d, s active
        // R = 3 * (4/9 + 1/9 + 1/9) = 3 * 6/9 = 2
        let r = r_ratio(2.0);
        assert!((r - 2.0).abs() < 1e-10, "R = {}", r);
    }

    #[test]
    fn test_r_ratio_above_charm() {
        // Above charm, below bottom: u, d, s, c active
        // R = 3 * (4/9 + 1/9 + 1/9 + 4/9) = 3 * 10/9 = 10/3
        let r = r_ratio(3.0);
        assert!((r - 10.0 / 3.0).abs() < 1e-10, "R = {}", r);
    }

    #[test]
    fn test_r_ratio_above_bottom() {
        // Above bottom, below top: u, d, s, c, b active
        // R = 3 * (4/9 + 1/9 + 1/9 + 4/9 + 1/9) = 3 * 11/9 = 11/3
        let r = r_ratio(10.0);
        assert!((r - 11.0 / 3.0).abs() < 1e-10, "R = {}", r);
    }
}
