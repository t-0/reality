use std::f64::consts::PI;

use crate::couplings::constants::GEV2_TO_NB;
use crate::kinematics::invariants::{kallen, Mandelstam};
use crate::matrix_elements::MatrixElement;

/// Compute the differential cross section d(sigma)/d(cos_theta) for a 2->2 process
/// at a given center-of-mass energy sqrt(s) and scattering angle cos(theta).
///
/// Returns the cross section in nanobarns.
///
/// Formula:
///   d(sigma)/d(cos_theta) = 1/(64*pi*s) * |p_f|/|p_i| * |M|^2_avg
///
/// For massless initial state (e+e-):
///   |p_i| = sqrt(s)/2
pub fn dsigma_dcostheta<M: MatrixElement>(
    me: &M,
    sqrt_s: f64,
    cos_theta: f64,
) -> f64 {
    let s = sqrt_s * sqrt_s;
    let (m3, m4) = me.final_masses();

    // Check threshold
    if sqrt_s <= m3 + m4 {
        return 0.0;
    }

    // Initial state momentum (massless approximation for e+e-)
    let p_i = sqrt_s / 2.0;

    // Final state momentum
    let p_f = kallen(s, m3 * m3, m4 * m4).max(0.0).sqrt() / (2.0 * sqrt_s);

    // Compute Mandelstam t from cos(theta)
    let mandelstam = Mandelstam::from_s_costheta(s, cos_theta, [0.0, 0.0], [m3, m4]);

    // Spin-averaged matrix element squared
    let m2 = me.amplitude_squared(s, mandelstam.t);

    // Differential cross section d(sigma)/d(cos_theta) (in GeV^-2)
    // = 2*pi * d(sigma)/d(Omega) = 2*pi * |M|^2 * |p_f| / (64*pi^2*s*|p_i|)
    //   = |M|^2 * |p_f| / (32*pi*s*|p_i|)
    let dsigma = 1.0 / (32.0 * PI * s) * (p_f / p_i) * m2;

    // Convert to nanobarns
    dsigma * GEV2_TO_NB
}

/// Compute d(sigma)/d(Omega) = d(sigma)/d(cos_theta) / (2*pi).
pub fn dsigma_domega<M: MatrixElement>(
    me: &M,
    sqrt_s: f64,
    cos_theta: f64,
) -> f64 {
    dsigma_dcostheta(me, sqrt_s, cos_theta) / (2.0 * PI)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::couplings::constants::ALPHA_EM;
    use crate::matrix_elements::ee_mumu::EpEmToMuPMuM;

    #[test]
    fn test_qed_angular_distribution() {
        // For e+e- -> mu+mu- QED, dsigma/dcostheta ~ (1 + cos^2 theta)
        let me = EpEmToMuPMuM::qed_only();
        let sqrt_s = 10.0;

        let forward = dsigma_dcostheta(&me, sqrt_s, 1.0);
        let perpendicular = dsigma_dcostheta(&me, sqrt_s, 0.0);

        // At cos_theta = 1: 1+1 = 2
        // At cos_theta = 0: 1+0 = 1
        // Ratio should be 2
        let ratio = forward / perpendicular;
        assert!(
            (ratio - 2.0).abs() < 0.01,
            "Angular ratio = {}, expected 2",
            ratio
        );
    }

    #[test]
    fn test_dsigma_positive() {
        let me = EpEmToMuPMuM::qed_only();
        let sqrt_s = 10.0;
        for cos_theta in [-0.9, -0.5, 0.0, 0.5, 0.9] {
            let ds = dsigma_dcostheta(&me, sqrt_s, cos_theta);
            assert!(ds > 0.0, "dsigma should be positive at cos_theta={}", cos_theta);
        }
    }
}
