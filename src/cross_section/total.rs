use crate::cross_section::integrator::gauss_legendre;
use crate::matrix_elements::MatrixElement;

use super::differential::dsigma_dcostheta;

/// Compute the total cross section for a 2->2 process by integrating
/// the differential cross section over cos(theta) from -1 to 1.
///
/// Returns the cross section in nanobarns.
///
/// Uses Gauss-Legendre quadrature with 64 points for high accuracy.
pub fn total_cross_section_2to2<M: MatrixElement>(me: &M, sqrt_s: f64) -> f64 {
    let (m3, m4) = me.final_masses();
    if sqrt_s <= m3 + m4 {
        return 0.0;
    }

    let n_points = 64;
    gauss_legendre(
        |cos_theta| dsigma_dcostheta(me, sqrt_s, cos_theta),
        -1.0,
        1.0,
        n_points,
    )
}

/// Compute the total cross section with a specified number of integration points.
pub fn total_cross_section_2to2_with_points<M: MatrixElement>(
    me: &M,
    sqrt_s: f64,
    n_points: usize,
) -> f64 {
    let (m3, m4) = me.final_masses();
    if sqrt_s <= m3 + m4 {
        return 0.0;
    }

    gauss_legendre(
        |cos_theta| dsigma_dcostheta(me, sqrt_s, cos_theta),
        -1.0,
        1.0,
        n_points,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::couplings::constants::*;
    use crate::matrix_elements::ee_mumu::EpEmToMuPMuM;
    use std::f64::consts::PI;

    #[test]
    fn test_total_xsec_qed_ee_mumu() {
        let me = EpEmToMuPMuM::qed_only();
        let sqrt_s = 10.0;
        let sigma = total_cross_section_2to2(&me, sqrt_s);

        // Analytic result: sigma = 4*pi*alpha^2/(3*s) in GeV^-2, then convert to nb
        let s = sqrt_s * sqrt_s;
        let sigma_expected = 4.0 * PI * ALPHA_EM * ALPHA_EM / (3.0 * s) * GEV2_TO_NB;

        assert!(
            (sigma - sigma_expected).abs() / sigma_expected < 0.01,
            "sigma = {} nb, expected = {} nb",
            sigma,
            sigma_expected
        );
    }

    #[test]
    fn test_z_pole_enhancement() {
        let me_ew = EpEmToMuPMuM::full_electroweak();
        let me_qed = EpEmToMuPMuM::qed_only();

        let sigma_ew = total_cross_section_2to2(&me_ew, M_Z);
        let sigma_qed = total_cross_section_2to2(&me_qed, M_Z);

        assert!(
            sigma_ew > 10.0 * sigma_qed,
            "Z pole should enhance cross section: EW = {} nb, QED = {} nb",
            sigma_ew,
            sigma_qed
        );
    }

    #[test]
    fn test_sigma_scales_as_1_over_s() {
        let me = EpEmToMuPMuM::qed_only();
        let sigma_10 = total_cross_section_2to2(&me, 10.0);
        let sigma_20 = total_cross_section_2to2(&me, 20.0);

        // sigma ~ 1/s, so sigma(20)/sigma(10) ~ (10/20)^2 = 0.25
        let ratio = sigma_20 / sigma_10;
        assert!(
            (ratio - 0.25).abs() < 0.01,
            "ratio = {}, expected ~0.25",
            ratio
        );
    }
}
