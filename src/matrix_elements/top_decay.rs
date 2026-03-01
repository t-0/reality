use std::f64::consts::PI;

use crate::couplings::constants::*;
use crate::particles::PdgId;

use super::amplitude::MatrixElement;

/// Top quark decay: t -> W+ b.
///
/// In the Standard Model, the top decays almost exclusively to Wb.
/// Partial width:
///   Gamma(t -> Wb) = (G_F * m_t^3) / (8 * pi * sqrt(2))
///                    * |V_tb|^2
///                    * (1 - m_W^2/m_t^2)^2 * (1 + 2*m_W^2/m_t^2)
/// (neglecting m_b)
pub struct TopDecay;

impl TopDecay {
    /// Analytic partial width for t -> W+ b (neglecting m_b).
    pub fn partial_width(&self) -> f64 {
        let mt = M_TOP;
        let mw = M_W;
        let vtb2 = 0.999118_f64.powi(2); // |V_tb|^2

        let r = mw * mw / (mt * mt); // r = m_W^2 / m_t^2

        G_FERMI * mt * mt * mt / (8.0 * PI * 2.0_f64.sqrt())
            * vtb2
            * (1.0 - r).powi(2)
            * (1.0 + 2.0 * r)
    }

    /// Width including finite m_b correction.
    pub fn partial_width_with_mb(&self) -> f64 {
        let mt = M_TOP;
        let mw = M_W;
        let mb = M_BOTTOM;
        let vtb2 = 0.999118_f64.powi(2);

        let x_w = mw * mw / (mt * mt);
        let x_b = mb * mb / (mt * mt);

        // Kallen function for the phase space
        let lambda = 1.0 + x_w * x_w + x_b * x_b - 2.0 * x_w - 2.0 * x_b - 2.0 * x_w * x_b;
        if lambda <= 0.0 {
            return 0.0;
        }

        let sqrt_lambda = lambda.sqrt();

        // Full width with m_b:
        // Gamma = GF * mt^3 / (8*pi*sqrt(2)) * |Vtb|^2 * sqrt(lambda) *
        //         [(1-xb)^2 + xw*(1+xb) - 2*xw^2]
        G_FERMI * mt * mt * mt / (8.0 * PI * 2.0_f64.sqrt())
            * vtb2
            * sqrt_lambda
            * ((1.0 - x_b).powi(2) + x_w * (1.0 + x_b) - 2.0 * x_w * x_w)
    }
}

impl MatrixElement for TopDecay {
    fn amplitude_squared(&self, s: f64, _t: f64) -> f64 {
        // s = m_t^2 for the decay
        // |M|^2 summed over final spins, averaged over initial spin (factor 1/2)
        // |M|^2 = g^2/(2*mw^2) * |V_tb|^2 * (2*p_t.p_b * p_W.p_t - m_t^2 * p_W.p_b)
        // In the rest frame with massless b:
        let mt2 = s;
        let mw2 = M_W * M_W;
        let vtb2 = 0.999118_f64.powi(2);
        let g2 = 4.0 * PI * ALPHA_EM / SIN2_THETA_W;

        // Simplified: |M|^2_avg = g^2 * |Vtb|^2 * mt2 * (1 - mw2/mt2)^2 * (1 + 2*mw2/mt2)
        g2 * vtb2 * mt2 * (1.0 - mw2 / mt2).powi(2) * (1.0 + 2.0 * mw2 / mt2)
    }

    fn initial_state(&self) -> Vec<PdgId> {
        vec![PdgId::TOP]
    }

    fn final_state(&self) -> Vec<PdgId> {
        vec![PdgId::W_BOSON, PdgId::BOTTOM]
    }

    fn initial_state_dof(&self) -> f64 {
        2.0 // spin-1/2 average: 2 spin states
    }

    fn final_masses(&self) -> (f64, f64) {
        (M_W, M_BOTTOM)
    }

    fn name(&self) -> &str {
        "t -> W+ b"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_top_width() {
        let decay = TopDecay;
        let width = decay.partial_width();
        // Expected: ~1.4 GeV
        assert!(
            (width - 1.4).abs() < 0.3,
            "Gamma(t -> Wb) = {} GeV",
            width
        );
    }

    #[test]
    fn test_top_width_with_mb() {
        let decay = TopDecay;
        let width = decay.partial_width_with_mb();
        // Should be close to the massless-b result
        let width_massless = decay.partial_width();
        let diff = (width - width_massless).abs() / width_massless;
        assert!(
            diff < 0.05,
            "m_b correction too large: {:.1}%",
            diff * 100.0
        );
    }

    #[test]
    fn test_top_lifetime() {
        let decay = TopDecay;
        let width = decay.partial_width();
        // Lifetime tau = hbar / Gamma
        // In natural units: tau = 1/Gamma (in GeV^-1)
        // Convert: 1 GeV^-1 = 6.582e-25 s
        let tau_gev = 1.0 / width;
        let tau_seconds = tau_gev * 6.582e-25;
        // Expected: ~5e-25 s
        assert!(
            tau_seconds > 1e-25 && tau_seconds < 1e-24,
            "tau_t = {} s",
            tau_seconds
        );
    }
}
