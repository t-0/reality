use std::f64::consts::PI;

use super::constants::{ALPHA_S_MZ, M_BOTTOM, M_CHARM, M_TOP, M_Z};

/// Number of active quark flavors at energy scale Q (in GeV).
pub fn n_active_flavors(q: f64) -> u32 {
    if q >= M_TOP {
        6
    } else if q >= M_BOTTOM {
        5
    } else if q >= M_CHARM {
        4
    } else {
        3
    }
}

/// One-loop QCD beta function coefficient: b0 = (33 - 2*n_f) / (12*pi).
fn beta0(n_f: u32) -> f64 {
    (33.0 - 2.0 * n_f as f64) / (12.0 * PI)
}

/// Running strong coupling constant alpha_s(Q) at one-loop.
///
/// Uses the reference value alpha_s(M_Z) = 0.1180 and runs to the given scale Q.
/// Crosses flavor thresholds at m_c, m_b, m_t.
pub fn alpha_s(q: f64) -> f64 {
    if q <= 0.0 {
        return 0.0;
    }

    // Start from alpha_s(M_Z)
    let mut alpha = ALPHA_S_MZ;
    let mut mu = M_Z;
    let target_nf = n_active_flavors(q);
    let start_nf = n_active_flavors(M_Z); // = 5

    // If we need to run to lower scales, step down through thresholds
    if target_nf < start_nf {
        let thresholds = [M_TOP, M_BOTTOM, M_CHARM];
        for &thresh in &thresholds {
            let nf = n_active_flavors(thresh + 0.001);
            if nf <= target_nf {
                break;
            }
            if mu > thresh && q < thresh {
                alpha = run_alpha_s(alpha, mu, thresh, nf);
                mu = thresh;
            }
        }
    }
    // If we need to run to higher scales, step up through thresholds
    else if target_nf > start_nf {
        if q >= M_TOP && mu < M_TOP {
            alpha = run_alpha_s(alpha, mu, M_TOP, 5);
            mu = M_TOP;
        }
    }

    run_alpha_s(alpha, mu, q, target_nf)
}

/// Run alpha_s from scale mu to scale q with n_f active flavors (one-loop).
/// alpha_s(q) = alpha_s(mu) / (1 + alpha_s(mu) * b0 * ln(q^2/mu^2))
fn run_alpha_s(alpha_mu: f64, mu: f64, q: f64, n_f: u32) -> f64 {
    let b0 = beta0(n_f);
    let log_ratio = (q * q / (mu * mu)).ln();
    let denom = 1.0 + alpha_mu * b0 * log_ratio;
    if denom <= 0.0 {
        // Landau pole — return a large value
        return 10.0;
    }
    alpha_mu / denom
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alpha_s_at_mz() {
        let a = alpha_s(M_Z);
        assert!((a - ALPHA_S_MZ).abs() < 1e-10);
    }

    #[test]
    fn test_alpha_s_increases_at_low_q() {
        // alpha_s should increase as Q decreases (asymptotic freedom in reverse)
        let a_high = alpha_s(100.0);
        let a_low = alpha_s(10.0);
        assert!(a_low > a_high, "alpha_s(10) = {} should be > alpha_s(100) = {}", a_low, a_high);
    }

    #[test]
    fn test_alpha_s_at_10_gev() {
        let a = alpha_s(10.0);
        // Expected ~ 0.18 at 10 GeV
        assert!(a > 0.15 && a < 0.25, "alpha_s(10 GeV) = {}", a);
    }

    #[test]
    fn test_alpha_s_decreases_at_high_q() {
        // Asymptotic freedom: alpha_s decreases at high Q
        let a_mz = alpha_s(M_Z);
        let a_1000 = alpha_s(1000.0);
        assert!(a_1000 < a_mz);
    }

    #[test]
    fn test_n_active_flavors() {
        assert_eq!(n_active_flavors(1.0), 3);
        assert_eq!(n_active_flavors(2.0), 4);
        assert_eq!(n_active_flavors(10.0), 5);
        assert_eq!(n_active_flavors(200.0), 6);
    }
}
