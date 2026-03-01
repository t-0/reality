
use num_complex::Complex64;

use super::vertex::VertexFactor;
use crate::couplings::constants::*;
use crate::dirac::DiracMatrix;

/// Charged current vertex: f-f'-W.
///
/// Feynman rule: -i * (g / 2*sqrt(2)) * gamma^mu * (1 - gamma5) * V_CKM
///
/// For leptons, V_CKM = 1 (no mixing).
pub fn charged_current_vertex(mu: usize, ckm_element: f64) -> VertexFactor {
    let g = g_weak();
    let coupling = -g / (2.0 * 2.0_f64.sqrt()) * ckm_element;

    // gamma^mu * (1 - gamma5) = gamma^mu * 2 * P_L
    let gamma_mu = DiracMatrix::gamma(mu);
    let pl = DiracMatrix::proj_left();
    let gamma_structure = gamma_mu * pl * 2.0;

    VertexFactor::FermionBoson {
        coupling: Complex64::new(coupling, 0.0),
        gamma_structure,
    }
}

/// Neutral current vertex: f-f-Z.
///
/// Feynman rule: -i * (g / cos(theta_W)) * gamma^mu * (g_V - g_A * gamma5) / 2
///
/// g_V = T3 - 2 * Q * sin^2(theta_W)
/// g_A = T3
pub fn neutral_current_vertex(mu: usize, t3: f64, charge: f64) -> VertexFactor {
    let g = g_weak();
    let gv = nc_vector_coupling(t3, charge);
    let ga = nc_axial_coupling(t3);

    let coupling = -g / (2.0 * COS_THETA_W);

    // gamma^mu * (g_V - g_A * gamma5)
    let gamma_mu = DiracMatrix::gamma(mu);
    let gamma5 = DiracMatrix::gamma5();
    let gamma_structure = gamma_mu * DiracMatrix::IDENTITY * gv - gamma_mu * gamma5 * ga;

    VertexFactor::FermionBoson {
        coupling: Complex64::new(coupling, 0.0),
        gamma_structure,
    }
}

/// Vector coupling g_V for a fermion with given T3 and charge.
pub fn gv_fermion(t3: f64, charge: f64) -> f64 {
    nc_vector_coupling(t3, charge)
}

/// Axial coupling g_A for a fermion with given T3.
pub fn ga_fermion(t3: f64) -> f64 {
    nc_axial_coupling(t3)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_electron_couplings() {
        // Electron: T3 = -1/2, Q = -1
        let gv = gv_fermion(-0.5, -1.0);
        let ga = ga_fermion(-0.5);
        // gV = -1/2 - 2*(-1)*0.23122 = -1/2 + 0.46244 = -0.03756
        assert!((gv - (-0.03756)).abs() < 0.001, "gV_e = {}", gv);
        assert!((ga - (-0.5)).abs() < 1e-10, "gA_e = {}", ga);
    }

    #[test]
    fn test_up_quark_couplings() {
        // Up quark: T3 = +1/2, Q = +2/3
        let gv = gv_fermion(0.5, 2.0 / 3.0);
        // gV = 1/2 - 2*(2/3)*0.23122 = 0.5 - 0.30829 = 0.19171
        assert!((gv - 0.19171).abs() < 0.001, "gV_u = {}", gv);
    }
}
