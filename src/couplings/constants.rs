use std::f64::consts::PI;

/// Fine structure constant at q^2 = 0.
pub const ALPHA_EM: f64 = 1.0 / 137.035999084;

/// Fine structure constant at the Z mass scale.
pub const ALPHA_EM_MZ: f64 = 1.0 / 127.952;

/// Strong coupling constant at M_Z.
pub const ALPHA_S_MZ: f64 = 0.1180;

/// Fermi coupling constant (GeV^-2).
pub const G_FERMI: f64 = 1.1663788e-5;

/// sin^2(theta_W) — Weinberg angle.
pub const SIN2_THETA_W: f64 = 0.23122;

/// cos^2(theta_W)
pub const COS2_THETA_W: f64 = 1.0 - SIN2_THETA_W;

/// sin(theta_W)
pub const SIN_THETA_W: f64 = 0.48076; // sqrt(0.23122)

/// cos(theta_W)
pub const COS_THETA_W: f64 = 0.87693; // sqrt(1 - 0.23122)

/// Electroweak vacuum expectation value (GeV).
/// v = 1 / sqrt(sqrt(2) * G_F) ~ 246.22 GeV
pub const V_EW: f64 = 246.22;

/// Weak coupling constant g_w = e / sin(theta_W)
pub fn g_weak() -> f64 {
    let e = (4.0 * PI * ALPHA_EM).sqrt();
    e / SIN_THETA_W
}

/// Weak neutral coupling g_z = e / (sin(theta_W) * cos(theta_W))
pub fn g_z() -> f64 {
    let e = (4.0 * PI * ALPHA_EM).sqrt();
    e / (SIN_THETA_W * COS_THETA_W)
}

/// Electric charge unit e = sqrt(4 * pi * alpha).
pub fn e_charge() -> f64 {
    (4.0 * PI * ALPHA_EM).sqrt()
}

/// Conversion factor: GeV^{-2} to picobarns.
pub const GEV2_TO_PB: f64 = 0.3894e9;

/// Conversion factor: GeV^{-2} to femtobarns.
pub const GEV2_TO_FB: f64 = 0.3894e12;

/// Conversion factor: GeV^{-2} to nanobarns.
pub const GEV2_TO_NB: f64 = 0.3894e6;

/// Z boson mass (GeV).
pub const M_Z: f64 = 91.1876;

/// W boson mass (GeV).
pub const M_W: f64 = 80.377;

/// Higgs boson mass (GeV).
pub const M_H: f64 = 125.20;

/// Top quark mass (GeV).
pub const M_TOP: f64 = 172.57;

/// Bottom quark mass (GeV).
pub const M_BOTTOM: f64 = 4.18;

/// Charm quark mass (GeV).
pub const M_CHARM: f64 = 1.27;

/// Tau lepton mass (GeV).
pub const M_TAU: f64 = 1.777;

/// Number of colors N_c.
pub const N_C: f64 = 3.0;

/// Color factor C_F = (N_c^2 - 1) / (2 * N_c) = 4/3.
pub const C_F: f64 = 4.0 / 3.0;

/// Color factor T_R = 1/2.
pub const T_R: f64 = 0.5;

/// Vector and axial couplings for neutral current: f-f-Z.
/// g_V^f = T3_f - 2 * Q_f * sin^2(theta_W)
/// g_A^f = T3_f
pub fn nc_vector_coupling(t3: f64, charge: f64) -> f64 {
    t3 - 2.0 * charge * SIN2_THETA_W
}

pub fn nc_axial_coupling(t3: f64) -> f64 {
    t3
}
