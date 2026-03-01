use std::f64::consts::PI;

use num_complex::Complex64;

use super::vertex::VertexFactor;
use crate::couplings::constants::ALPHA_S_MZ;
use crate::dirac::DiracMatrix;

/// QCD quark-quark-gluon vertex.
///
/// Feynman rule: -i * g_s * gamma^mu * T^a
/// For color-summed calculations, the color factor is handled separately.
pub fn qqg_vertex(mu: usize, alpha_s: f64) -> VertexFactor {
    let g_s = (4.0 * PI * alpha_s).sqrt();
    VertexFactor::FermionBoson {
        coupling: Complex64::new(-g_s, 0.0),
        gamma_structure: DiracMatrix::gamma(mu),
    }
}

/// Strong coupling g_s at the Z mass scale.
pub fn g_strong_mz() -> f64 {
    (4.0 * PI * ALPHA_S_MZ).sqrt()
}

/// Triple gluon vertex coupling.
/// V = g_s * f^{abc} * [g^{mu nu}(k1-k2)^rho + cyclic permutations]
pub fn triple_gluon_coupling(alpha_s: f64) -> VertexFactor {
    let g_s = (4.0 * PI * alpha_s).sqrt();
    VertexFactor::ThreeBoson {
        coupling: Complex64::new(g_s, 0.0),
    }
}

/// Quartic gluon vertex coupling.
/// V = -i * g_s^2 * [color and Lorentz structure]
pub fn quartic_gluon_coupling(alpha_s: f64) -> VertexFactor {
    let g_s = (4.0 * PI * alpha_s).sqrt();
    VertexFactor::FourBoson {
        coupling: Complex64::new(-g_s * g_s, 0.0),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gs_mz() {
        let gs = g_strong_mz();
        // g_s(M_Z) = sqrt(4 * pi * 0.118) ~ 1.218
        assert!((gs - 1.218).abs() < 0.01, "g_s = {}", gs);
    }
}
