use num_complex::Complex64;

use crate::dirac::DiracMatrix;

/// Represents a Feynman rule vertex factor.
#[derive(Debug, Clone)]
pub enum VertexFactor {
    /// Fermion-fermion-boson vertex: coupling * gamma_structure.
    /// The full vertex is: i * coupling * gamma_structure
    FermionBoson {
        coupling: Complex64,
        gamma_structure: DiracMatrix,
    },
    /// Three-boson vertex (e.g., triple gauge coupling).
    /// The Lorentz structure depends on the external momenta and is handled
    /// at the matrix-element level.
    ThreeBoson { coupling: Complex64 },
    /// Four-boson vertex.
    FourBoson { coupling: Complex64 },
    /// Scalar-fermion-fermion (Yukawa) vertex.
    Yukawa { coupling: Complex64 },
}
