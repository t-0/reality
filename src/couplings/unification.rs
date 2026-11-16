use std::f64::consts::PI;

use super::constants::{ALPHA_EM_MZ, ALPHA_S_MZ, COS2_THETA_W, M_Z, SIN2_THETA_W};

/// One-loop beta coefficients (b1, b2, b3) for the Standard Model,
/// with U(1)_Y in GUT normalization (g1^2 = 5/3 g'^2).
pub const B_SM: [f64; 3] = [41.0 / 10.0, -19.0 / 6.0, -7.0];

/// One-loop beta coefficients (b1, b2, b3) for the MSSM.
pub const B_MSSM: [f64; 3] = [33.0 / 5.0, 1.0, -3.0];

/// Upper limit of the scale scan (GeV), just below the Planck scale.
pub const Q_MAX: f64 = 1.0e19;

/// Inverse gauge couplings (1/alpha_1, 1/alpha_2, 1/alpha_3) at M_Z.
///
/// alpha_1 = (5/3) alpha_em / cos^2(theta_W)   (SU(5)/SO(10) normalization)
/// alpha_2 = alpha_em / sin^2(theta_W)
/// alpha_3 = alpha_s
pub fn inverse_couplings_mz() -> [f64; 3] {
    [
        0.6 * COS2_THETA_W / ALPHA_EM_MZ,
        SIN2_THETA_W / ALPHA_EM_MZ,
        1.0 / ALPHA_S_MZ,
    ]
}

/// Run inverse couplings at one loop from scale mu to scale q:
/// 1/alpha_i(q) = 1/alpha_i(mu) - b_i / (2*pi) * ln(q/mu)
fn run_inverse(inv_mu: [f64; 3], mu: f64, q: f64, b: [f64; 3]) -> [f64; 3] {
    let l = (q / mu).ln() / (2.0 * PI);
    [inv_mu[0] - b[0] * l, inv_mu[1] - b[1] * l, inv_mu[2] - b[2] * l]
}

/// Inverse gauge couplings at scale q (GeV), running from M_Z.
///
/// With `susy_scale = None` the pure Standard Model is used. With
/// `Some(m_susy)` the SM runs up to m_susy and the MSSM above it
/// (a single, sharp superpartner threshold; values below M_Z are clamped).
pub fn inverse_couplings(q: f64, susy_scale: Option<f64>) -> [f64; 3] {
    let inv_mz = inverse_couplings_mz();
    match susy_scale {
        Some(ms) if q > ms.max(M_Z) => {
            let ms = ms.max(M_Z);
            let inv_ms = run_inverse(inv_mz, M_Z, ms, B_SM);
            run_inverse(inv_ms, ms, q, B_MSSM)
        }
        _ => run_inverse(inv_mz, M_Z, q, B_SM),
    }
}

/// The scale where two couplings meet.
#[derive(Debug, Clone, Copy)]
pub struct Crossing {
    /// Indices (0-based) of the two couplings: 0 = U(1), 1 = SU(2), 2 = SU(3).
    pub pair: (usize, usize),
    /// Scale in GeV.
    pub scale: f64,
    /// Common 1/alpha at the crossing.
    pub inverse_alpha: f64,
}

/// Find where couplings i and j cross between M_Z and Q_MAX, by bisection in ln Q.
pub fn crossing(i: usize, j: usize, susy_scale: Option<f64>) -> Option<Crossing> {
    let diff = |ln_q: f64| {
        let a = inverse_couplings(ln_q.exp(), susy_scale);
        a[i] - a[j]
    };
    let (mut lo, mut hi) = (M_Z.ln(), Q_MAX.ln());
    if diff(lo).signum() == diff(hi).signum() {
        return None;
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if diff(lo).signum() == diff(mid).signum() {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let scale = (0.5 * (lo + hi)).exp();
    Some(Crossing {
        pair: (i, j),
        scale,
        inverse_alpha: inverse_couplings(scale, susy_scale)[i],
    })
}

/// All three pairwise crossings: (1,2), (2,3), (1,3).
pub fn crossings(susy_scale: Option<f64>) -> [Option<Crossing>; 3] {
    [
        crossing(0, 1, susy_scale),
        crossing(1, 2, susy_scale),
        crossing(0, 2, susy_scale),
    ]
}

/// The scale where the three couplings come closest, measured by the
/// spread max(1/alpha_i) - min(1/alpha_i). Returns (scale, spread).
pub fn closest_approach(susy_scale: Option<f64>) -> (f64, f64) {
    let (lo, hi) = (M_Z.ln(), Q_MAX.ln());
    let n = 4000;
    let mut best = (M_Z, f64::INFINITY);
    for k in 0..=n {
        let q = (lo + (hi - lo) * k as f64 / n as f64).exp();
        let a = inverse_couplings(q, susy_scale);
        let spread = a.iter().cloned().fold(f64::MIN, f64::max)
            - a.iter().cloned().fold(f64::MAX, f64::min);
        if spread < best.1 {
            best = (q, spread);
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_inputs_at_mz() {
        let a = inverse_couplings_mz();
        assert!((a[0] - 59.0).abs() < 0.5, "1/alpha_1 = {}", a[0]);
        assert!((a[1] - 29.6).abs() < 0.5, "1/alpha_2 = {}", a[1]);
        assert!((a[2] - 8.47).abs() < 0.05, "1/alpha_3 = {}", a[2]);
    }

    #[test]
    fn test_sm_does_not_unify() {
        // SM: alpha_1 = alpha_2 near 1e13 GeV, alpha_2 = alpha_3 near 1e17 GeV.
        let c12 = crossing(0, 1, None).unwrap();
        let c23 = crossing(1, 2, None).unwrap();
        assert!(c12.scale > 1e12 && c12.scale < 1e14, "{:e}", c12.scale);
        assert!(c23.scale > 1e16 && c23.scale < 1e18, "{:e}", c23.scale);
        let (_, spread) = closest_approach(None);
        assert!(spread > 2.0, "SM spread = {}", spread);
    }

    #[test]
    fn test_mssm_unifies() {
        // MSSM from ~1 TeV: all three meet near 2e16 GeV with 1/alpha ~ 25.
        let c12 = crossing(0, 1, Some(1000.0)).unwrap();
        let c23 = crossing(1, 2, Some(1000.0)).unwrap();
        assert!(c12.scale > 5e15 && c12.scale < 5e16, "{:e}", c12.scale);
        assert!((c12.scale / c23.scale).ln().abs() < 1.0);
        assert!((c12.inverse_alpha - 25.0).abs() < 3.0);
        let (_, spread) = closest_approach(Some(1000.0));
        assert!(spread < 1.0, "MSSM spread = {}", spread);
    }

    #[test]
    fn test_sin2w_at_gut_scale_is_three_eighths() {
        // When alpha_1 = alpha_2, sin^2(theta_W) = 3/8 in GUT normalization.
        let c = crossing(0, 1, None).unwrap();
        let a = inverse_couplings(c.scale, None);
        // sin^2 = g'^2 / (g^2 + g'^2) with alpha' = (3/5) alpha_1.
        let (alpha1, alpha2) = (1.0 / a[0], 1.0 / a[1]);
        let sin2 = 0.6 * alpha1 / (alpha2 + 0.6 * alpha1);
        assert!((sin2 - 0.375).abs() < 1e-6, "sin2w = {}", sin2);
    }
}
