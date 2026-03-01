use num_complex::Complex64;

use super::vertex::VertexFactor;
use crate::couplings::constants::V_EW;

/// Yukawa vertex: f-fbar-H.
///
/// Feynman rule: -i * (m_f / v) where v = 246.22 GeV.
pub fn yukawa_vertex(fermion_mass: f64) -> VertexFactor {
    let coupling = -fermion_mass / V_EW;
    VertexFactor::Yukawa {
        coupling: Complex64::new(coupling, 0.0),
    }
}

/// Yukawa coupling for a fermion: y_f = sqrt(2) * m_f / v
pub fn yukawa_coupling(fermion_mass: f64) -> f64 {
    2.0_f64.sqrt() * fermion_mass / V_EW
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::couplings::constants::M_TOP;

    #[test]
    fn test_top_yukawa() {
        let y_t = yukawa_coupling(M_TOP);
        // y_t = sqrt(2) * 172.57 / 246.22 ~ 0.991
        assert!((y_t - 0.991).abs() < 0.01, "y_t = {}", y_t);
    }

    #[test]
    fn test_yukawa_vertex() {
        if let VertexFactor::Yukawa { coupling } = yukawa_vertex(M_TOP) {
            assert!((coupling.re - (-M_TOP / V_EW)).abs() < 1e-10);
        } else {
            panic!("Wrong vertex type");
        }
    }
}
