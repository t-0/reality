use std::f64::consts::PI;

use crate::couplings::constants::*;
use crate::particles::PdgId;

use super::amplitude::MatrixElement;

/// H -> f fbar tree-level decay.
///
/// Partial width: Gamma(H -> f fbar) = N_c * m_f^2 * m_H / (8 * pi * v^2) * beta^3
/// where beta = sqrt(1 - 4*m_f^2/m_H^2) and N_c is the color factor.
///
/// The matrix element is:
///   |M|^2 = N_c * (m_f / v)^2 * (m_H^2 - 4*m_f^2)
pub struct HiggsToFermions {
    pub fermion_pdg: PdgId,
    pub fermion_mass: f64,
    pub n_colors: f64,
}

impl HiggsToFermions {
    pub fn new(pdg: PdgId, mass: f64, n_c: f64) -> Self {
        Self {
            fermion_pdg: pdg,
            fermion_mass: mass,
            n_colors: n_c,
        }
    }

    pub fn to_bottom() -> Self {
        Self::new(PdgId::BOTTOM, M_BOTTOM, N_C)
    }

    pub fn to_charm() -> Self {
        Self::new(PdgId::CHARM, M_CHARM, N_C)
    }

    pub fn to_tau() -> Self {
        Self::new(PdgId::TAU, M_TAU, 1.0) // lepton: N_c = 1
    }

    /// Analytic partial width in GeV.
    pub fn partial_width(&self) -> f64 {
        let mh = M_H;
        let mf = self.fermion_mass;

        if mh < 2.0 * mf {
            return 0.0;
        }

        let beta = (1.0 - 4.0 * mf * mf / (mh * mh)).sqrt();
        self.n_colors * mf * mf * mh / (8.0 * PI * V_EW * V_EW) * beta * beta * beta
    }
}

impl MatrixElement for HiggsToFermions {
    fn amplitude_squared(&self, s: f64, _t: f64) -> f64 {
        // s = m_H^2 for decay
        let mf = self.fermion_mass;
        // |M|^2_avg over initial Higgs spin (trivial, spin-0) and summed over final spins
        // |M|^2 = N_c * (m_f/v)^2 * 2 * (s - 4*m_f^2)
        self.n_colors * (mf / V_EW).powi(2) * 2.0 * (s - 4.0 * mf * mf)
    }

    fn initial_state(&self) -> Vec<PdgId> {
        vec![PdgId::HIGGS]
    }

    fn final_state(&self) -> Vec<PdgId> {
        vec![self.fermion_pdg, self.fermion_pdg.antiparticle()]
    }

    fn initial_state_dof(&self) -> f64 {
        1.0 // scalar
    }

    fn final_masses(&self) -> (f64, f64) {
        (self.fermion_mass, self.fermion_mass)
    }

    fn name(&self) -> &str {
        "H -> f fbar"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_higgs_to_bb_width() {
        let h_bb = HiggsToFermions::to_bottom();
        let width = h_bb.partial_width();
        // Tree-level with pole mass m_b = 4.18 GeV gives ~4.3 MeV.
        // The PDG value ~2.4 MeV uses running mass m_b(m_H) ~ 2.8 GeV with QCD corrections.
        assert!(
            width > 0.003 && width < 0.006,
            "Gamma(H->bb) = {} GeV (tree-level, pole mass)",
            width
        );
    }

    #[test]
    fn test_higgs_to_tau_width() {
        let h_tt = HiggsToFermions::to_tau();
        let width = h_tt.partial_width();
        // Expected: ~0.26 MeV = 0.00026 GeV
        assert!(
            width > 0.0001 && width < 0.001,
            "Gamma(H->tau tau) = {} GeV",
            width
        );
    }

    #[test]
    fn test_higgs_to_cc_width() {
        let h_cc = HiggsToFermions::to_charm();
        let width = h_cc.partial_width();
        // Expected: ~0.12 MeV
        assert!(width > 0.0 && width < 0.001, "Gamma(H->cc) = {} GeV", width);
    }
}
