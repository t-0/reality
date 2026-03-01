use num_complex::Complex64;

use super::gamma::DiracMatrix;
use crate::kinematics::FourMomentum;

/// Compute the trace of a product of DiracMatrix values.
pub fn trace_product(matrices: &[DiracMatrix]) -> Complex64 {
    if matrices.is_empty() {
        return DiracMatrix::IDENTITY.trace();
    }
    let mut product = matrices[0];
    for m in &matrices[1..] {
        product = product * *m;
    }
    product.trace()
}

/// Compute Tr[gamma^mu gamma^nu] = 4 g^{mu nu}
pub fn trace_two_gamma(mu: usize, nu: usize) -> f64 {
    if mu == nu {
        let metric = [1.0, -1.0, -1.0, -1.0];
        4.0 * metric[mu]
    } else {
        0.0
    }
}

/// Compute the spin-summed amplitude-squared for a fermion line,
/// using trace technology:
///   sum_spins |ubar(p3) * Gamma * u(p1)|^2
///     = Tr[(p3_slash + m3) * Gamma * (p1_slash + m1) * Gamma_bar]
/// where Gamma_bar = gamma0 * Gamma^dagger * gamma0.
pub fn fermion_trace(
    p_out: &FourMomentum,
    m_out: f64,
    vertex_gamma: &DiracMatrix,
    p_in: &FourMomentum,
    m_in: f64,
) -> Complex64 {
    let sum_out = DiracMatrix::slash(p_out) + DiracMatrix::IDENTITY * m_out;
    let sum_in = DiracMatrix::slash(p_in) + DiracMatrix::IDENTITY * m_in;
    let gamma_bar = vertex_gamma.adjoint();
    (sum_out * *vertex_gamma * sum_in * gamma_bar).trace()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trace_product_identity() {
        assert!((trace_product(&[]) - Complex64::new(4.0, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_trace_product_single_gamma() {
        // Tr[gamma^mu] = 0
        for mu in 0..4 {
            let tr = trace_product(&[DiracMatrix::gamma(mu)]);
            assert!(tr.norm() < 1e-10);
        }
    }

    #[test]
    fn test_trace_product_two_gammas() {
        let metric = [1.0, -1.0, -1.0, -1.0];
        for mu in 0..4 {
            for nu in 0..4 {
                let tr = trace_product(&[DiracMatrix::gamma(mu), DiracMatrix::gamma(nu)]);
                let expected = if mu == nu { 4.0 * metric[mu] } else { 0.0 };
                assert!(
                    (tr - Complex64::new(expected, 0.0)).norm() < 1e-10,
                    "Tr[g{} g{}] = {:?}, expected {}",
                    mu, nu, tr, expected
                );
            }
        }
    }

    #[test]
    fn test_trace_odd_gammas_zero() {
        // Tr[gamma^mu gamma^nu gamma^rho] = 0
        let tr = trace_product(&[
            DiracMatrix::gamma(0),
            DiracMatrix::gamma(1),
            DiracMatrix::gamma(2),
        ]);
        assert!(tr.norm() < 1e-10);
    }

    #[test]
    fn test_trace_four_gammas() {
        // Tr[gamma^mu gamma^nu gamma^rho gamma^sigma]
        //   = 4 * (g^{mn} g^{rs} - g^{mr} g^{ns} + g^{ms} g^{nr})
        let g = |mu: usize, nu: usize| -> f64 {
            if mu == nu {
                [1.0, -1.0, -1.0, -1.0][mu]
            } else {
                0.0
            }
        };

        for mu in 0..4 {
            for nu in 0..4 {
                for rho in 0..4 {
                    for sigma in 0..4 {
                        let tr = trace_product(&[
                            DiracMatrix::gamma(mu),
                            DiracMatrix::gamma(nu),
                            DiracMatrix::gamma(rho),
                            DiracMatrix::gamma(sigma),
                        ]);
                        let expected = 4.0
                            * (g(mu, nu) * g(rho, sigma) - g(mu, rho) * g(nu, sigma)
                                + g(mu, sigma) * g(nu, rho));
                        assert!(
                            (tr - Complex64::new(expected, 0.0)).norm() < 1e-10,
                            "4-gamma trace failed: mu={}, nu={}, rho={}, sigma={}",
                            mu,
                            nu,
                            rho,
                            sigma
                        );
                    }
                }
            }
        }
    }
}
