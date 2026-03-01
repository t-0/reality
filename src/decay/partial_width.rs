use std::f64::consts::PI;

use crate::couplings::constants::*;

/// Partial width for W -> l nu_l (one lepton flavor).
/// Gamma(W -> l nu) = G_F * M_W^3 / (6 * pi * sqrt(2))
pub fn w_to_lepton_nu() -> f64 {
    G_FERMI * M_W * M_W * M_W / (6.0 * PI * 2.0_f64.sqrt())
}

/// Partial width for W -> q_i q_j_bar (one quark channel).
/// Gamma(W -> qi qj) = N_c * G_F * M_W^3 / (6 * pi * sqrt(2)) * |V_ij|^2
/// Neglecting quark masses for light quarks.
pub fn w_to_quarks(v_ij_squared: f64) -> f64 {
    N_C * G_FERMI * M_W * M_W * M_W / (6.0 * PI * 2.0_f64.sqrt()) * v_ij_squared
}

/// Partial width for Z -> f fbar.
/// Gamma(Z -> f fbar) = N_c * G_F * M_Z^3 / (6 * pi * sqrt(2)) * (g_V^2 + g_A^2)
pub fn z_to_ffbar(t3: f64, charge: f64, n_c: f64) -> f64 {
    let gv = nc_vector_coupling(t3, charge);
    let ga = nc_axial_coupling(t3);
    n_c * G_FERMI * M_Z * M_Z * M_Z / (6.0 * PI * 2.0_f64.sqrt()) * (gv * gv + ga * ga)
}

/// Partial width for H -> f fbar.
/// Gamma(H -> f fbar) = N_c * m_f^2 * M_H / (8 * pi * v^2) * beta^3
pub fn h_to_ffbar(fermion_mass: f64, n_c: f64) -> f64 {
    let mh = M_H;
    if mh < 2.0 * fermion_mass {
        return 0.0;
    }
    let beta = (1.0 - 4.0 * fermion_mass * fermion_mass / (mh * mh)).sqrt();
    n_c * fermion_mass * fermion_mass * mh / (8.0 * PI * V_EW * V_EW) * beta.powi(3)
}

/// Partial width for t -> W+ b.
pub fn top_to_wb() -> f64 {
    let mt = M_TOP;
    let mw = M_W;
    let vtb2 = VTB_SQUARED;
    let r = mw * mw / (mt * mt);
    G_FERMI * mt * mt * mt / (8.0 * PI * 2.0_f64.sqrt()) * vtb2 * (1.0 - r).powi(2) * (1.0 + 2.0 * r)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_w_to_enu() {
        let w = w_to_lepton_nu();
        // Expected: ~227 MeV = 0.227 GeV
        assert!(
            (w - 0.227).abs() < 0.02,
            "Gamma(W->e nu) = {} GeV",
            w
        );
    }

    #[test]
    fn test_z_to_ee() {
        let w = z_to_ffbar(-0.5, -1.0, 1.0);
        // Expected: ~84 MeV = 0.084 GeV
        assert!(
            (w - 0.084).abs() < 0.005,
            "Gamma(Z->ee) = {} GeV",
            w
        );
    }
}
