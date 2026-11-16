use std::f64::consts::PI;

use num_complex::Complex64;

use crate::couplings::constants::*;
use crate::dirac::gamma::DiracMatrix;
use crate::dirac::spinor::{DiracSpinor, SpinState};
use crate::kinematics::polarization::{conjugate, dot_cc, dot_rc, polarization_vectors};
use crate::kinematics::FourMomentum;
use crate::particles::PdgId;

use super::amplitude::MatrixElement;

/// e+e- -> W+W- via s-channel photon, s-channel Z, and t-channel neutrino exchange.
///
/// Uses exact numerical helicity amplitude evaluation with the full Dirac
/// spinor and polarization vector infrastructure. All three diagrams are
/// summed coherently at the amplitude level for each spin/polarization
/// configuration.
pub struct EpEmToWpWm;

impl MatrixElement for EpEmToWpWm {
    fn amplitude_squared(&self, s: f64, t: f64) -> f64 {
        let mw = M_W;
        let mw2 = mw * mw;
        let sqrt_s = s.sqrt();

        // Check threshold
        if sqrt_s < 2.0 * mw {
            return 0.0;
        }

        let e_beam = sqrt_s / 2.0;
        let beta = (1.0 - 4.0 * mw2 / s).sqrt();

        if beta < 1e-12 {
            return 0.0;
        }

        let p_w = e_beam * beta;

        // Reconstruct cos_theta from s and t
        // t = mw² - s/2 + (s·β/2)·cosθ  →  cosθ = (2t + s - 2mw²) / (s·β)
        let cos_theta = ((2.0 * t + s - 2.0 * mw2) / (s * beta)).clamp(-1.0, 1.0);
        let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();

        // CM frame 4-momenta (massless electrons)
        let p1 = FourMomentum::new(e_beam, 0.0, 0.0, e_beam); // e⁻
        let p2 = FourMomentum::new(e_beam, 0.0, 0.0, -e_beam); // e⁺
        let k1 = FourMomentum::new(e_beam, p_w * sin_theta, 0.0, p_w * cos_theta); // W⁻
        let k2 = FourMomentum::new(e_beam, -p_w * sin_theta, 0.0, -p_w * cos_theta); // W⁺

        // Neutrino propagator momentum
        let p_nu = p1 - k1;
        let t_val = p_nu.invariant_mass_sq(); // = (p₁ - k₁)² < 0

        // Coupling constants using α(M_Z)
        let e2 = 4.0 * PI * ALPHA_EM_MZ;
        let sw2 = SIN2_THETA_W;

        // Z propagator: 1/(s - M_Z² + i M_Z Γ_Z)
        let mz2 = M_Z * M_Z;
        let d_z = Complex64::new(1.0, 0.0) / Complex64::new(s - mz2, M_Z * GAMMA_Z);

        // Electron neutral current couplings for Z
        let gv_e = nc_vector_coupling(-0.5, -1.0);
        let ga_e = nc_axial_coupling(-0.5);

        // Precompute Z coupling matrix: (g_V·I - g_A·γ⁵)
        let z_coupling_base =
            DiracMatrix::IDENTITY.scale_real(gv_e) - DiracMatrix::gamma5().scale_real(ga_e);

        // Polarization vectors for W⁻ and W⁺
        let eps1_all = polarization_vectors(&k1, mw);
        let eps2_all = polarization_vectors(&k2, mw);

        // Precompute momentum combinations for the triple gauge vertex
        // Using all-incoming convention: W⁻(-k₁), W⁺(-k₂), V(q=k₁+k₂)
        // Vertex momenta: (k₂-k₁), (-k₁-2k₂), (2k₁+k₂)
        let k2_minus_k1 = k2 - k1;
        let neg_k1_minus_2k2 = FourMomentum::new(
            -k1.e - 2.0 * k2.e,
            -k1.px - 2.0 * k2.px,
            -k1.py - 2.0 * k2.py,
            -k1.pz - 2.0 * k2.pz,
        );
        let two_k1_plus_k2 = FourMomentum::new(
            2.0 * k1.e + k2.e,
            2.0 * k1.px + k2.px,
            2.0 * k1.py + k2.py,
            2.0 * k1.pz + k2.pz,
        );

        let metric: [f64; 4] = [1.0, -1.0, -1.0, -1.0];
        let spins = [SpinState::Up, SpinState::Down];

        let mut sum_m2 = 0.0;

        for &s1 in &spins {
            for &s2 in &spins {
                // Dirac spinors for massless electrons
                let u = DiracSpinor::u(&p1, 0.0, s1);
                let v = DiracSpinor::v(&p2, 0.0, s2);
                let vbar = v.bar();

                // Compute lepton currents J^μ for photon and Z
                let mut j_gamma = [Complex64::new(0.0, 0.0); 4];
                let mut j_z = [Complex64::new(0.0, 0.0); 4];

                for mu in 0..4 {
                    let gamma_mu = DiracMatrix::gamma(mu);
                    j_gamma[mu] = vbar.sandwich(&gamma_mu, &u);

                    let z_vertex = gamma_mu.multiply(&z_coupling_base);
                    j_z[mu] = vbar.sandwich(&z_vertex, &u);
                }

                for lam1 in 0..3 {
                    for lam2 in 0..3 {
                        // Complex conjugate polarization vectors (outgoing particles)
                        let eps1c = conjugate(&eps1_all[lam1]);
                        let eps2c = conjugate(&eps2_all[lam2]);

                        // Triple gauge vertex contraction T^μ (same structure for γ and Z):
                        // Using all-incoming convention:
                        // T^μ = (ε₁*·ε₂*)(k₂-k₁)^μ + [(-k₁-2k₂)·ε₁*]ε₂*^μ
                        //        + [(2k₁+k₂)·ε₂*]ε₁*^μ
                        let eps12_dot = dot_cc(&eps1c, &eps2c);
                        let coeff_eps2 = dot_rc(&neg_k1_minus_2k2, &eps1c);
                        let coeff_eps1 = dot_rc(&two_k1_plus_k2, &eps2c);

                        let k12_arr = k2_minus_k1.as_array();

                        let mut jt_gamma = Complex64::new(0.0, 0.0);
                        let mut jt_z = Complex64::new(0.0, 0.0);

                        for mu in 0..4 {
                            let t_mu = eps12_dot * Complex64::new(k12_arr[mu], 0.0)
                                + coeff_eps2 * eps2c[mu]
                                + coeff_eps1 * eps1c[mu];

                            let m = metric[mu];
                            jt_gamma += j_gamma[mu] * t_mu * m;
                            jt_z += j_z[mu] * t_mu * m;
                        }

                        // s-channel photon: M_γ = (Q_e e²/s) × J_γ · T with Q_e = -1.
                        // The sign is what lets the photon cancel the Q-dependent part
                        // of the Z coupling at high energy.
                        let m_gamma = Complex64::new(-e2 / s, 0.0) * jt_gamma;

                        // s-channel Z: M_Z = [e²/(2sin²θ_W)] × D_Z × J_Z · T
                        let m_z = Complex64::new(e2 / (2.0 * sw2), 0.0) * d_z * jt_z;

                        // t-channel neutrino: M_ν = [e²/(2sin²θ_W·t)] × v̄ ε̸₂* P_L p̸_ν ε̸₁* P_L u
                        let chain = DiracMatrix::slash_complex(&eps2c)
                            .multiply(&DiracMatrix::proj_left())
                            .multiply(&DiracMatrix::slash(&p_nu))
                            .multiply(&DiracMatrix::slash_complex(&eps1c))
                            .multiply(&DiracMatrix::proj_left());
                        let nu_spinor = vbar.sandwich(&chain, &u);
                        let m_nu =
                            Complex64::new(e2 / (2.0 * sw2 * t_val), 0.0) * nu_spinor;

                        let m_total = m_gamma + m_z + m_nu;
                        sum_m2 += m_total.norm_sqr();
                    }
                }
            }
        }

        // Spin average: divide by 4 (2 electron spins × 2 positron spins)
        sum_m2 / 4.0
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cross_section::differential::dsigma_dcostheta;
    use crate::cross_section::total::total_cross_section_2to2;

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
        let mw2 = M_W * M_W;
        let t = mw2 - s / 2.0; // cos_theta = 0
        let m2 = me.amplitude_squared(s, t);
        assert!(m2 > 0.0, "|M|^2 should be positive, got {}", m2);
    }

    #[test]
    fn test_total_cross_section_200gev() {
        // Tree-level SM prediction at √s = 200 GeV is approximately 17 pb
        let me = EpEmToWpWm;
        let sigma_nb = total_cross_section_2to2(&me, 200.0);
        let sigma_pb = sigma_nb * 1000.0;

        assert!(
            sigma_pb > 14.0 && sigma_pb < 20.0,
            "σ(200 GeV) = {:.2} pb, expected ~17 pb",
            sigma_pb
        );
    }

    #[test]
    fn test_gauge_cancellation_high_energy() {
        // Without gauge cancellation, σ would grow with s for longitudinal W's.
        // The SM prediction should decrease roughly as 1/s at high energy.
        let me = EpEmToWpWm;

        let sigma_200 = total_cross_section_2to2(&me, 200.0);
        let sigma_500 = total_cross_section_2to2(&me, 500.0);
        let sigma_1000 = total_cross_section_2to2(&me, 1000.0);

        assert!(
            sigma_500 > 0.0,
            "σ(500 GeV) should be positive: {}",
            sigma_500
        );
        assert!(
            sigma_1000 > 0.0,
            "σ(1000 GeV) should be positive: {}",
            sigma_1000
        );
        assert!(
            sigma_500 < sigma_200,
            "σ should decrease with energy: σ(500)={} > σ(200)={}",
            sigma_500,
            sigma_200
        );
        assert!(
            sigma_1000 < sigma_500,
            "σ should decrease with energy: σ(1000)={} > σ(500)={}",
            sigma_1000,
            sigma_500
        );
    }

    #[test]
    fn test_threshold_behavior() {
        // Near threshold the t-channel neutrino dominates and sigma rises linearly
        // in beta (S-wave): sigma -> pi alpha^2 beta / (sin^4(theta_W) s).
        let me = EpEmToWpWm;
        let leading = |sqrt_s: f64| {
            let s = sqrt_s * sqrt_s;
            let beta = (1.0 - 4.0 * M_W * M_W / s).sqrt();
            PI * ALPHA_EM_MZ * ALPHA_EM_MZ * beta / (SIN2_THETA_W * SIN2_THETA_W * s) * GEV2_TO_NB
        };
        for de in [0.05, 0.5] {
            let sqrt_s = 2.0 * M_W + de;
            let sigma = total_cross_section_2to2(&me, sqrt_s);
            let expected = leading(sqrt_s);
            assert!(
                (sigma / expected - 1.0).abs() < 0.03,
                "sqrt(s) = 2 M_W + {} GeV: sigma = {} nb, leading-order threshold formula = {} nb",
                de,
                sigma,
                expected
            );
        }
        assert_eq!(total_cross_section_2to2(&me, 2.0 * M_W - 0.1), 0.0);
    }

    #[test]
    fn debug_diagram_contributions() {
        let sqrt_s = 500.0_f64;
        let s = sqrt_s * sqrt_s;
        let mw = M_W;
        let mw2 = mw * mw;
        let e_beam = sqrt_s / 2.0;
        let beta = (1.0 - 4.0 * mw2 / s).sqrt();
        let p_w = e_beam * beta;
        let cos_theta = 0.0_f64;
        let sin_theta = 1.0_f64;
        let e2 = 4.0 * PI * ALPHA_EM_MZ;
        let sw2 = SIN2_THETA_W;
        let mz2 = M_Z * M_Z;
        let d_z = Complex64::new(1.0, 0.0) / Complex64::new(s - mz2, M_Z * GAMMA_Z);
        let gv_e = nc_vector_coupling(-0.5, -1.0);
        let ga_e = nc_axial_coupling(-0.5);

        let p1 = FourMomentum::new(e_beam, 0.0, 0.0, e_beam);
        let p2 = FourMomentum::new(e_beam, 0.0, 0.0, -e_beam);
        let k1 = FourMomentum::new(e_beam, p_w * sin_theta, 0.0, p_w * cos_theta);
        let k2 = FourMomentum::new(e_beam, -p_w * sin_theta, 0.0, -p_w * cos_theta);
        let p_nu = p1 - k1;
        let t_val = p_nu.invariant_mass_sq();

        eprintln!("sqrt_s={}, beta={:.4}, p_w={:.2}, t={:.2}", sqrt_s, beta, p_w, t_val);
        eprintln!("e2={:.6}, e2/(2sw2)={:.6}, 1/s={:.2e}, d_z={:?}", e2, e2/(2.0*sw2), 1.0/s, d_z);

        let eps1_all = polarization_vectors(&k1, mw);
        let eps2_all = polarization_vectors(&k2, mw);
        let z_coupling_base = DiracMatrix::IDENTITY.scale_real(gv_e) - DiracMatrix::gamma5().scale_real(ga_e);

        let k2_minus_k1 = k2 - k1;
        let neg_k1_minus_2k2 = FourMomentum::new(-k1.e-2.0*k2.e, -k1.px-2.0*k2.px, -k1.py-2.0*k2.py, -k1.pz-2.0*k2.pz);
        let two_k1_plus_k2 = FourMomentum::new(2.0*k1.e+k2.e, 2.0*k1.px+k2.px, 2.0*k1.py+k2.py, 2.0*k1.pz+k2.pz);
        let metric_diag: [f64; 4] = [1.0, -1.0, -1.0, -1.0];

        let mut sum_gamma2 = 0.0_f64;
        let mut sum_z2 = 0.0_f64;
        let mut sum_nu2 = 0.0_f64;
        let mut sum_total2 = 0.0_f64;

        for &s1 in &[SpinState::Up, SpinState::Down] {
            for &s2 in &[SpinState::Up, SpinState::Down] {
                let u = DiracSpinor::u(&p1, 0.0, s1);
                let v = DiracSpinor::v(&p2, 0.0, s2);
                let vbar = v.bar();

                let mut j_gamma = [Complex64::new(0.0,0.0); 4];
                let mut j_z = [Complex64::new(0.0,0.0); 4];
                for mu in 0..4 {
                    let gmu = DiracMatrix::gamma(mu);
                    j_gamma[mu] = vbar.sandwich(&gmu, &u);
                    j_z[mu] = vbar.sandwich(&gmu.multiply(&z_coupling_base), &u);
                }

                for lam1 in 0..3 {
                    for lam2 in 0..3 {
                        let eps1c = conjugate(&eps1_all[lam1]);
                        let eps2c = conjugate(&eps2_all[lam2]);

                        let eps12_dot = dot_cc(&eps1c, &eps2c);
                        let coeff_eps2 = dot_rc(&neg_k1_minus_2k2, &eps1c);
                        let coeff_eps1 = dot_rc(&two_k1_plus_k2, &eps2c);
                        let k12_arr = k2_minus_k1.as_array();

                        let mut jt_gamma = Complex64::new(0.0, 0.0);
                        let mut jt_z = Complex64::new(0.0, 0.0);
                        for mu in 0..4 {
                            let t_mu = eps12_dot * Complex64::new(k12_arr[mu], 0.0) + coeff_eps2 * eps2c[mu] + coeff_eps1 * eps1c[mu];
                            jt_gamma += j_gamma[mu] * t_mu * metric_diag[mu];
                            jt_z += j_z[mu] * t_mu * metric_diag[mu];
                        }

                        let m_gamma = Complex64::new(-e2 / s, 0.0) * jt_gamma;
                        let m_z = Complex64::new(e2 / (2.0 * sw2), 0.0) * d_z * jt_z;

                        let chain = DiracMatrix::slash_complex(&eps2c)
                            .multiply(&DiracMatrix::proj_left())
                            .multiply(&DiracMatrix::slash(&p_nu))
                            .multiply(&DiracMatrix::slash_complex(&eps1c))
                            .multiply(&DiracMatrix::proj_left());
                        let nu_spinor = vbar.sandwich(&chain, &u);
                        let m_nu = Complex64::new(e2 / (2.0 * sw2 * t_val), 0.0) * nu_spinor;

                        let m_total = m_gamma + m_z + m_nu;

                        sum_gamma2 += m_gamma.norm_sqr();
                        sum_z2 += m_z.norm_sqr();
                        sum_nu2 += m_nu.norm_sqr();
                        sum_total2 += m_total.norm_sqr();

                        // Print longitudinal-longitudinal (lam=1 is ε(0))
                        if lam1 == 1 && lam2 == 1 {
                            eprintln!("spin=({:?},{:?}) LL: γ={:.4e} Z={:.4e} ν={:.4e} tot={:.4e}",
                                s1, s2, m_gamma, m_z, m_nu, m_total);
                        }
                    }
                }
            }
        }

        eprintln!("\nΣ|M|²/4: γ={:.6e}, Z={:.6e}, ν={:.6e}, total={:.6e}",
            sum_gamma2/4.0, sum_z2/4.0, sum_nu2/4.0, sum_total2/4.0);

        // Convert to cross section
        let p_i = e_beam;
        let p_f = p_w;
        let avg_m2 = sum_total2 / 4.0;
        let sigma_gev2 = 1.0/(32.0 * PI * s) * (p_f/p_i) * avg_m2 * 2.0; // approx integral over cosθ
        let sigma_pb = sigma_gev2 * GEV2_TO_NB * 1000.0;
        eprintln!("Approx σ ~ {:.2} pb (very rough)", sigma_pb);
    }

    #[test]
    fn test_angular_distribution_forward_backward() {
        // t-channel neutrino exchange enhances forward scattering
        // (W⁻ goes in the e⁻ direction)
        let me = EpEmToWpWm;
        let sqrt_s = 200.0;

        let forward = dsigma_dcostheta(&me, sqrt_s, 0.5);
        let backward = dsigma_dcostheta(&me, sqrt_s, -0.5);

        assert!(
            forward > backward,
            "Forward ({}) should exceed backward ({}) at 200 GeV",
            forward,
            backward
        );
    }
}
