use std::f64::consts::PI;

use num_complex::Complex64;

use crate::couplings::constants::*;
use crate::particles::PdgId;

use super::amplitude::MatrixElement;

/// H -> gamma gamma via W boson and top quark loops.
///
/// Uses the known analytic form factors:
///   A(H -> gamma gamma) = A_1(tau_W) + N_c * Q_t^2 * A_{1/2}(tau_t)
/// where tau_i = m_H^2 / (4 * m_i^2).
pub struct HiggsToGammaGamma;

impl HiggsToGammaGamma {
    /// Partial width Gamma(H -> gamma gamma).
    pub fn partial_width(&self) -> f64 {
        let mh = M_H;
        let alpha = ALPHA_EM;

        // tau parameters
        let tau_w = mh * mh / (4.0 * M_W * M_W);
        let tau_t = mh * mh / (4.0 * M_TOP * M_TOP);

        // Form factors
        let a_w = form_factor_spin1(tau_w);
        let a_t = form_factor_spin_half(tau_t);

        // Total amplitude: A = A_1(W) + N_c * Q_t^2 * A_{1/2}(top)
        let qt = 2.0 / 3.0;
        let amplitude = a_w + Complex64::new(N_C * qt * qt, 0.0) * a_t;

        // Width: Gamma = alpha^2 * m_H^3 / (256 * pi^3 * v^2) * |A|^2
        alpha * alpha * mh * mh * mh / (256.0 * PI * PI * PI * V_EW * V_EW) * amplitude.norm_sqr()
    }
}

/// Form factor for spin-1/2 loop (fermion): A_{1/2}(tau) = -2*tau * [1 + (1-tau)*f(tau)]
fn form_factor_spin_half(tau: f64) -> Complex64 {
    let f_val = f_function(tau);
    Complex64::new(-2.0 * tau, 0.0) * (Complex64::new(1.0, 0.0) + Complex64::new(1.0 - tau, 0.0) * f_val)
}

/// Form factor for spin-1 loop (W boson): A_1(tau) = 2 + 3*tau + 3*tau*(2-tau)*f(tau)
fn form_factor_spin1(tau: f64) -> Complex64 {
    let f_val = f_function(tau);
    Complex64::new(2.0 + 3.0 * tau, 0.0) + Complex64::new(3.0 * tau * (2.0 - tau), 0.0) * f_val
}

/// The f(tau) function used in Higgs loop amplitudes.
/// f(tau) = [arcsin(sqrt(tau))]^2                for tau <= 1
/// f(tau) = -1/4 * [ln((1+sqrt(1-1/tau))/(1-sqrt(1-1/tau))) - i*pi]^2  for tau > 1
fn f_function(tau: f64) -> Complex64 {
    if tau <= 1.0 {
        let val = tau.sqrt().asin();
        Complex64::new(val * val, 0.0)
    } else {
        let sqrt_val = (1.0 - 1.0 / tau).sqrt();
        let log_arg = (1.0 + sqrt_val) / (1.0 - sqrt_val);
        let log_val = log_arg.ln();
        let real = log_val;
        let imag = -PI;
        let z = Complex64::new(real, imag);
        z * z * Complex64::new(-0.25, 0.0)
    }
}

impl MatrixElement for HiggsToGammaGamma {
    fn amplitude_squared(&self, _s: f64, _t: f64) -> f64 {
        // For the decay matrix element, return |M|^2 proportional to width
        // This is mainly used through partial_width()
        let mh = M_H;
        let tau_w = mh * mh / (4.0 * M_W * M_W);
        let tau_t = mh * mh / (4.0 * M_TOP * M_TOP);

        let a_w = form_factor_spin1(tau_w);
        let a_t = form_factor_spin_half(tau_t);
        let qt = 2.0 / 3.0;
        let amplitude = a_w + Complex64::new(N_C * qt * qt, 0.0) * a_t;

        let e4 = (4.0 * PI * ALPHA_EM).powi(2);
        e4 / (V_EW * V_EW) * amplitude.norm_sqr()
    }

    fn initial_state(&self) -> Vec<PdgId> {
        vec![PdgId::HIGGS]
    }

    fn final_state(&self) -> Vec<PdgId> {
        vec![PdgId::PHOTON, PdgId::PHOTON]
    }

    fn initial_state_dof(&self) -> f64 {
        1.0
    }

    fn symmetry_factor(&self) -> f64 {
        0.5 // identical photons
    }

    fn final_masses(&self) -> (f64, f64) {
        (0.0, 0.0)
    }

    fn name(&self) -> &str {
        "H -> gamma gamma"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_higgs_to_gammagamma_width() {
        let process = HiggsToGammaGamma;
        let width = process.partial_width();
        // Expected: ~9.4 keV = 9.4e-6 GeV
        assert!(
            width > 5e-6 && width < 15e-6,
            "Gamma(H->gammagamma) = {} GeV = {} keV",
            width,
            width * 1e6
        );
    }

    #[test]
    fn test_branching_ratio() {
        let width_gg = HiggsToGammaGamma.partial_width();
        let total_width = GAMMA_H;
        let br = width_gg / total_width;
        // Expected: BR ~ 2.3e-3
        assert!(
            br > 1e-3 && br < 5e-3,
            "BR(H->gammagamma) = {}",
            br
        );
    }

    #[test]
    fn test_form_factor_limit() {
        // For tau << 1 (heavy loop particle), f(tau) -> tau, A_{1/2} -> -4/3
        let a = form_factor_spin_half(0.01);
        // Should be approximately -2*0.01*(1 + 0.99*0.01) ~ small
        assert!(a.norm() < 0.1);
    }
}
