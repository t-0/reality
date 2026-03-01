use std::f64::consts::PI;

use crate::couplings::constants::*;
use crate::particles::PdgId;

use super::amplitude::MatrixElement;

/// e+e- -> W+W- via s-channel photon, s-channel Z, and t-channel neutrino exchange.
///
/// This is the most complex tree-level process. The three diagrams must be
/// summed coherently at the amplitude level.
///
/// Uses the high-energy approximation for the squared amplitude.
pub struct EpEmToWpWm;

impl MatrixElement for EpEmToWpWm {
    fn amplitude_squared(&self, s: f64, t: f64) -> f64 {
        let sqrt_s = s.sqrt();
        let mw = M_W;
        let mw2 = mw * mw;

        // Check threshold
        if sqrt_s < 2.0 * mw {
            return 0.0;
        }

        let u = 2.0 * mw2 - s - t;

        let e2 = 4.0 * PI * ALPHA_EM;
        let sw2 = SIN2_THETA_W;
        let _cw2 = COS2_THETA_W;

        // Propagators
        let prop_gamma = 1.0 / s; // photon s-channel
        let mz2 = M_Z * M_Z;
        let gamma_z = GAMMA_Z;
        let denom_z = (s - mz2) * (s - mz2) + (M_Z * gamma_z) * (M_Z * gamma_z);
        let prop_z_re = (s - mz2) / denom_z;
        let prop_z_im = -M_Z * gamma_z / denom_z;

        // t-channel neutrino propagator: 1/t (massless neutrino)
        let prop_nu = if t.abs() > 1e-30 { 1.0 / t } else { 0.0 };

        // CM momentum of W
        let beta_w = (1.0 - 4.0 * mw2 / s).max(0.0).sqrt();
        let _p_w = sqrt_s / 2.0 * beta_w;

        // Use the standard result for unpolarized cross section:
        // dsigma/dOmega = alpha^2 / (4s) * beta_w * F(s, cos_theta)
        //
        // The full expression involves a lengthy combination of diagrams.
        // Here we use the compact analytic form.
        //
        // F(s,t) from the three-diagram interference:

        // Diagram amplitudes (squared and interference):
        // A_gamma: s-channel photon
        // A_Z: s-channel Z
        // A_nu: t-channel neutrino

        // For the gamma-gamma term:
        let gg = e2 * e2 * prop_gamma * prop_gamma;
        // For the Z-Z term:
        let zz = e2 * e2 / (sw2 * sw2) * (prop_z_re * prop_z_re + prop_z_im * prop_z_im);
        // For the nu-nu term:
        let nn = e2 * e2 / (4.0 * sw2 * sw2) * prop_nu * prop_nu;
        // For gamma-Z interference (real part):
        let gz = e2 * e2 / sw2 * prop_gamma * prop_z_re;
        // For gamma-nu interference:
        let gn = e2 * e2 / (2.0 * sw2) * prop_gamma * prop_nu;
        // For Z-nu interference:
        let zn = e2 * e2 / (2.0 * sw2 * sw2) * prop_nu * prop_z_re;

        // Kinematic factors for each term
        // These come from contracting the polarization tensors and fermion traces
        // Using t_hat = t - mw2, u_hat = u - mw2
        let _t_hat = t - mw2;
        let _u_hat = u - mw2;

        // Simplified form using Mandelstam variables:
        // The full result for the spin-averaged amplitude squared:
        // See e.g. Peskin & Schroeder Problem 15.5, or
        // M. Bohm, A. Denner, H. Joos "Gauge Theories" eq. (8.3.40)

        // s-channel (gamma + Z combined) contribution:
        let f_s = gg * self.s_channel_factor(s, t, mw2)
            + zz * self.s_channel_z_factor(s, t, mw2, sw2)
            + gz * 2.0 * self.s_channel_gz_factor(s, t, mw2, sw2);

        // t-channel neutrino contribution:
        let f_t = nn * self.t_channel_factor(s, t, mw2);

        // s-t interference:
        let f_st = gn * self.st_gamma_nu_factor(s, t, mw2)
            + zn * self.st_z_nu_factor(s, t, mw2, sw2);

        f_s + f_t + 2.0 * f_st
    }

    fn initial_state(&self) -> Vec<PdgId> {
        vec![PdgId::POSITRON, PdgId::ELECTRON]
    }

    fn final_state(&self) -> Vec<PdgId> {
        vec![PdgId::W_BOSON, PdgId::W_MINUS]
    }

    fn initial_state_dof(&self) -> f64 {
        4.0
    }

    fn final_masses(&self) -> (f64, f64) {
        (M_W, M_W)
    }

    fn name(&self) -> &str {
        "e+e- -> W+W-"
    }
}

impl EpEmToWpWm {
    // s-channel photon tensor contraction
    fn s_channel_factor(&self, s: f64, t: f64, mw2: f64) -> f64 {
        let u = 2.0 * mw2 - s - t;
        // Photon couples to W charge. The result from contracting
        // the triple gauge vertex with fermion current:
        // ~ (t*u - mw^4) * (2 + (s - 2mw2)^2 / (4 mw^4)) + correction terms
        //
        // Use the standard Denner et al. result:
        // For unpolarized photon exchange:
        let tu = t * u;
        let s2 = s * s;

        // Transverse part
        let term1 = tu / mw2 - s2;
        // Longitudinal corrections
        let term2 = s2 + 8.0 * mw2 * s - 8.0 * mw2 * mw2;

        -(term1 + term2 / 4.0 - s * s / (4.0 * mw2))
    }

    fn s_channel_z_factor(&self, s: f64, t: f64, mw2: f64, _sw2: f64) -> f64 {
        // Z boson has additional (gV_e, gA_e) couplings and modified WWZ vertex
        let gv_e = nc_vector_coupling(-0.5, -1.0);
        let ga_e = nc_axial_coupling(-0.5);
        let v2_plus_a2 = gv_e * gv_e + ga_e * ga_e;

        v2_plus_a2 * self.s_channel_factor(s, t, mw2)
    }

    fn s_channel_gz_factor(&self, s: f64, t: f64, mw2: f64, _sw2: f64) -> f64 {
        let gv_e = nc_vector_coupling(-0.5, -1.0);
        gv_e * self.s_channel_factor(s, t, mw2)
    }

    fn t_channel_factor(&self, s: f64, t: f64, mw2: f64) -> f64 {
        let _u = 2.0 * mw2 - s - t;
        // t-channel neutrino exchange gives a left-handed structure
        // |M_nu|^2 ~ (t - mw2)^2 * (s/mw2 + ...) terms
        let t_hat = t - mw2;

        // Simplified expression for the t-channel squared amplitude
        4.0 * t_hat * t_hat / mw2 + 4.0 * t_hat + s
    }

    fn st_gamma_nu_factor(&self, s: f64, t: f64, mw2: f64) -> f64 {
        let t_hat = t - mw2;
        // gamma-nu interference
        -(2.0 * t_hat + s) * (1.0 + t_hat / mw2)
    }

    fn st_z_nu_factor(&self, s: f64, t: f64, mw2: f64, _sw2: f64) -> f64 {
        let gv_e = nc_vector_coupling(-0.5, -1.0);
        let ga_e = nc_axial_coupling(-0.5);
        // Z-nu interference involves both gV and gA
        (gv_e + ga_e) * self.st_gamma_nu_factor(s, t, mw2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_below_threshold() {
        let me = EpEmToWpWm;
        let s = (150.0_f64).powi(2); // 150 GeV < 2*M_W
        let t = -s / 2.0;
        assert_eq!(me.amplitude_squared(s, t), 0.0);
    }

    #[test]
    fn test_above_threshold_positive() {
        let me = EpEmToWpWm;
        let sqrt_s = 200.0;
        let s = sqrt_s * sqrt_s;
        // cos_theta = 0
        let mw2 = M_W * M_W;
        let p_w = (s / 4.0 - mw2).max(0.0).sqrt();
        let e_w = sqrt_s / 2.0;
        let t = mw2 - sqrt_s / 2.0 * (e_w - 0.0); // cos_theta=0, need proper formula
        // t = mw^2 - s/2 + sqrt(s) * p_w * cos_theta
        // For cos_theta = 0: t = mw^2 - s/2
        let t = mw2 - s / 2.0;

        let m2 = me.amplitude_squared(s, t);
        assert!(m2 > 0.0, "|M|^2 should be positive, got {}", m2);
    }
}
