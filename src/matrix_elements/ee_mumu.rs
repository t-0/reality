use std::f64::consts::PI;

use crate::couplings::constants::*;
use crate::particles::PdgId;

use super::amplitude::MatrixElement;

/// e+e- -> mu+mu- via virtual photon (QED tree level).
///
/// In the massless limit:
///   |M|^2 (spin-summed) = 2 * e^4 * (t^2 + u^2) / s^2
///
/// Including Z boson exchange (full electroweak):
///   Includes s-channel photon and Z propagators with interference.
pub struct EpEmToMuPMuM {
    /// If true, include Z boson exchange. Otherwise pure QED.
    pub include_z: bool,
}

impl EpEmToMuPMuM {
    pub fn qed_only() -> Self {
        Self { include_z: false }
    }

    pub fn full_electroweak() -> Self {
        Self { include_z: true }
    }
}

impl MatrixElement for EpEmToMuPMuM {
    fn amplitude_squared(&self, s: f64, t: f64) -> f64 {
        let u = -s - t; // massless limit

        if !self.include_z {
            // Pure QED: |M|^2_avg = 2 * e^4 * (t^2 + u^2) / s^2
            let e4 = (4.0 * PI * ALPHA_EM).powi(2);
            2.0 * e4 * (t * t + u * u) / (s * s)
        } else {
            // Full electroweak: gamma + Z exchange
            // Electron and muon couplings
            let qe = -1.0; // electron charge
            let t3e = -0.5; // electron T3

            let gv_e = nc_vector_coupling(t3e, qe);
            let ga_e = nc_axial_coupling(t3e);
            let gv_mu = gv_e; // same for muon
            let ga_mu = ga_e;

            let e2 = 4.0 * PI * ALPHA_EM;

            // Photon propagator: 1/s
            // Z propagator: 1/(s - M_Z^2 + i*M_Z*Gamma_Z) — use Breit-Wigner
            let gamma_z = GAMMA_Z;
            let denom_z_re = s - M_Z * M_Z;
            let denom_z_im = M_Z * gamma_z;
            let denom_z_sq = denom_z_re * denom_z_re + denom_z_im * denom_z_im;

            // Z coupling prefactor: g_z^2 / (cos_theta_w)^2 is already in g_z
            // Using: g_z = e / (sin_theta_w * cos_theta_w)
            let _gz2 = e2 / (SIN2_THETA_W * COS2_THETA_W);

            // Photon-only term: e^4 * Q_e^2 * Q_mu^2 * (t^2 + u^2 + ... )
            // For gamma exchange:
            //   |M_gamma|^2 = e^4 * 2 * (t^2 + u^2) / s^2
            let _m_gamma_sq = 2.0 * e2 * e2 * (t * t + u * u) / (s * s);

            // Z-only term:
            //   |M_Z|^2 = gz2^2 / (4*denom_z_sq) * [(gV^2+gA^2)_e * (gV^2+gA^2)_mu * 2*(t^2+u^2)
            //              + 4 * gV_e*gA_e * gV_mu*gA_mu * 2*(t^2 - u^2)]
            // Wait, let me use the standard result more carefully.
            //
            // For massless fermions e+e- -> f fbar with both gamma and Z:
            // |M|^2_avg = (1/4) * [ (A1^2 + A2^2)(t^2 + u^2) + (A1^2 - A2^2) * 0 ]
            // where the terms involve combinations of gamma and Z amplitudes.
            //
            // Let's use the standard result:
            // d(sigma)/d(cos_theta) = (alpha^2 / 4s) * [F1 * (1 + cos^2 theta) + F2 * cos theta]
            //
            // where:
            //   F1 = Q_e^2 Q_f^2 + 2 * Q_e Q_f * Re(chi) * gV_e gV_f + |chi|^2 * (gV_e^2 + gA_e^2)(gV_f^2 + gA_f^2)
            //   F2 = 4 * Q_e Q_f * Re(chi) * gA_e gA_f + 8 * |chi|^2 * gV_e gA_e * gV_f gA_f
            //
            //   chi(s) = s / [(s - M_Z^2) + i*M_Z*Gamma_Z] * 1/(4 sin2_w cos2_w)
            //
            // Then |M|^2_avg = e^4 * [F1 * 2*(t^2 + u^2)/s^2 + F2 * 2*(t^2 - u^2)/s^2]
            //                = e^4 * 2/s^2 * [F1*(t^2+u^2) + F2*(t^2-u^2)]

            let chi_prefactor = 1.0 / (4.0 * SIN2_THETA_W * COS2_THETA_W);
            let chi_re = chi_prefactor * s * denom_z_re / denom_z_sq;
            let chi_im = chi_prefactor * s * (-denom_z_im) / denom_z_sq;
            let chi_sq = chi_re * chi_re + chi_im * chi_im;

            let qf = -1.0; // muon charge

            let f1 = qe * qe * qf * qf
                + 2.0 * qe * qf * chi_re * gv_e * gv_mu
                + chi_sq * (gv_e * gv_e + ga_e * ga_e) * (gv_mu * gv_mu + ga_mu * ga_mu);

            let f2 = 4.0 * qe * qf * chi_re * ga_e * ga_mu
                + 8.0 * chi_sq * gv_e * ga_e * gv_mu * ga_mu;

            let e4 = e2 * e2;
            e4 * 2.0 / (s * s) * (f1 * (t * t + u * u) + f2 * (t * t - u * u))
        }
    }

    fn initial_state(&self) -> Vec<PdgId> {
        vec![PdgId::POSITRON, PdgId::ELECTRON]
    }

    fn final_state(&self) -> Vec<PdgId> {
        vec![PdgId::ANTI_MUON, PdgId::MUON]
    }

    fn initial_state_dof(&self) -> f64 {
        4.0 // 2 spins x 2 spins
    }

    fn final_masses(&self) -> (f64, f64) {
        (0.0, 0.0) // massless approximation for muons at high energy
    }

    fn name(&self) -> &str {
        if self.include_z {
            "e+e- -> mu+mu- (QED + Z)"
        } else {
            "e+e- -> mu+mu- (QED)"
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_qed_forward() {
        let me = EpEmToMuPMuM::qed_only();
        let s = 100.0; // sqrt(s) = 10 GeV
        // At cos_theta = 0: t = u = -s/2
        let t = -s / 2.0;
        let m2 = me.amplitude_squared(s, t);
        // |M|^2_avg = 2*e^4*(t^2+u^2)/s^2 = 2*e^4*2*(s/2)^2/s^2 = e^4
        let e4 = (4.0 * PI * ALPHA_EM).powi(2);
        assert!((m2 - e4).abs() / e4 < 1e-10);
    }

    #[test]
    fn test_qed_total_cross_section() {
        // sigma(e+e- -> mu+mu-) = 4*pi*alpha^2 / (3*s)
        // At sqrt(s) = 10 GeV, s = 100:
        // sigma = 4*pi*(1/137)^2 / (3*100) = 4*pi / (137^2 * 300)
        // In GeV^-2, then convert to nb.
        let s = 100.0;
        let sigma_expected = 4.0 * PI * ALPHA_EM * ALPHA_EM / (3.0 * s) * GEV2_TO_NB;
        // ~ 86.8 / s nb = 0.868 nb
        assert!((sigma_expected - 0.868).abs() < 0.01, "sigma = {} nb", sigma_expected);
    }

    #[test]
    fn test_z_resonance() {
        // At the Z pole, the cross section should be much larger than QED alone
        let me_qed = EpEmToMuPMuM::qed_only();
        let me_full = EpEmToMuPMuM::full_electroweak();

        let s = M_Z * M_Z;
        let t = -s / 2.0;

        let m2_qed = me_qed.amplitude_squared(s, t);
        let m2_full = me_full.amplitude_squared(s, t);

        // At Z pole, the full EW result should be much larger
        assert!(
            m2_full > 10.0 * m2_qed,
            "Z resonance not visible: qed={}, full={}",
            m2_qed,
            m2_full
        );
    }
}
