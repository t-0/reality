use std::f64::consts::PI;

use num_complex::Complex64;

use super::vertex::VertexFactor;
use crate::couplings::constants::ALPHA_EM;
use crate::dirac::DiracMatrix;

/// QED fermion-fermion-photon vertex.
///
/// Feynman rule: -i * e * Q_f * gamma^mu
///
/// The vertex factor returned is for a specific Lorentz index mu.
pub fn qed_vertex(charge: f64, mu: usize) -> VertexFactor {
    let e = (4.0 * PI * ALPHA_EM).sqrt();
    VertexFactor::FermionBoson {
        coupling: Complex64::new(-e * charge, 0.0),
        gamma_structure: DiracMatrix::gamma(mu),
    }
}

/// QED coupling constant e.
pub fn qed_coupling() -> f64 {
    (4.0 * PI * ALPHA_EM).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_electron_vertex_coupling() {
        let e = qed_coupling();
        // e ~ 0.3028
        assert!((e - 0.3028).abs() < 0.001, "e = {}", e);
    }

    #[test]
    fn test_qed_vertex_electron() {
        if let VertexFactor::FermionBoson { coupling, .. } = qed_vertex(-1.0, 0) {
            // For electron (Q = -1): coupling = -e * (-1) = +e
            let e = qed_coupling();
            assert!((coupling.re - e).abs() < 1e-10);
        } else {
            panic!("Wrong vertex type");
        }
    }
}
