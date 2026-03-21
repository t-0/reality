use num_complex::Complex64;

use crate::kinematics::FourMomentum;

/// Compute the three polarization vectors (λ = +1, 0, -1) for a massive spin-1
/// particle with the given 4-momentum and mass.
///
/// The momentum must lie in the xz-plane (py = 0).
///
/// Returns [ε(+1), ε(0), ε(-1)] as arrays of Complex64 contravariant components.
///
/// Conventions:
///   ε(±1) = ∓(1/√2)(0, n_z, ±i, -n_x)   (helicity basis)
///   ε(0)  = (|p|, E n_x, 0, E n_z) / M   (longitudinal)
///
/// where n̂ = (n_x, 0, n_z) is the unit vector along the 3-momentum.
pub fn polarization_vectors(k: &FourMomentum, mass: f64) -> [[Complex64; 4]; 3] {
    let p = k.three_momentum_mag();
    let e = k.e;

    let inv_sqrt2 = std::f64::consts::FRAC_1_SQRT_2;
    let zero = Complex64::new(0.0, 0.0);

    if p < 1e-10 {
        // Particle at rest: use z-axis as default direction
        // ε(+1) = -(1/√2)(0, 1, i, 0) → (0, -1/√2, -i/√2, 0)
        // ε(-1) = +(1/√2)(0, 1, -i, 0) → (0, 1/√2, -i/√2, 0)
        return [
            [
                zero,
                Complex64::new(-inv_sqrt2, 0.0),
                Complex64::new(0.0, -inv_sqrt2),
                zero,
            ],
            [zero, zero, zero, Complex64::new(1.0, 0.0)],
            [
                zero,
                Complex64::new(inv_sqrt2, 0.0),
                Complex64::new(0.0, -inv_sqrt2),
                zero,
            ],
        ];
    }

    let nx = k.px / p;
    let nz = k.pz / p;

    // ε(+1) = -(1/√2)(0, n_z, i, -n_x)
    let eps_plus = [
        zero,
        Complex64::new(-nz * inv_sqrt2, 0.0),
        Complex64::new(0.0, -inv_sqrt2),
        Complex64::new(nx * inv_sqrt2, 0.0),
    ];

    // ε(0) = (|p|, E·n_x, 0, E·n_z) / M
    let eps_zero = [
        Complex64::new(p / mass, 0.0),
        Complex64::new(e * nx / mass, 0.0),
        zero,
        Complex64::new(e * nz / mass, 0.0),
    ];

    // ε(-1) = +(1/√2)(0, n_z, -i, -n_x)  →  (0, n_z/√2, -i/√2, -n_x/√2)
    let eps_minus = [
        zero,
        Complex64::new(nz * inv_sqrt2, 0.0),
        Complex64::new(0.0, -inv_sqrt2),
        Complex64::new(-nx * inv_sqrt2, 0.0),
    ];

    [eps_plus, eps_zero, eps_minus]
}

/// Complex conjugate of a polarization 4-vector.
pub fn conjugate(eps: &[Complex64; 4]) -> [Complex64; 4] {
    [eps[0].conj(), eps[1].conj(), eps[2].conj(), eps[3].conj()]
}

/// Minkowski dot product of two complex 4-vectors:
/// a·b = a⁰b⁰ - a¹b¹ - a²b² - a³b³
pub fn dot_cc(a: &[Complex64; 4], b: &[Complex64; 4]) -> Complex64 {
    a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]
}

/// Minkowski dot product of a real 4-momentum with a complex 4-vector.
pub fn dot_rc(k: &FourMomentum, eps: &[Complex64; 4]) -> Complex64 {
    Complex64::new(k.e, 0.0) * eps[0]
        - Complex64::new(k.px, 0.0) * eps[1]
        - Complex64::new(k.py, 0.0) * eps[2]
        - Complex64::new(k.pz, 0.0) * eps[3]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polarization_orthogonality() {
        // W boson with some momentum
        let mw: f64 = 80.377;
        let e: f64 = 100.0;
        let p = (e * e - mw * mw).sqrt();
        let cos_theta: f64 = 0.6;
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
        let k = FourMomentum::new(e, p * sin_theta, 0.0, p * cos_theta);

        let eps = polarization_vectors(&k, mw);

        // Each polarization should be transverse to k: ε·k = 0
        for lam in 0..3 {
            let ek = dot_rc(&k, &eps[lam]);
            assert!(
                ek.norm() < 1e-10,
                "ε(λ={})·k = {:?}, expected 0",
                lam,
                ek
            );
        }

        // Normalization: ε(λ)·ε*(λ) = -1
        for lam in 0..3 {
            let ec = conjugate(&eps[lam]);
            let norm = dot_cc(&eps[lam], &ec);
            assert!(
                (norm.re + 1.0).abs() < 1e-10 && norm.im.abs() < 1e-10,
                "ε(λ={})·ε*(λ) = {:?}, expected -1",
                lam,
                norm
            );
        }

        // Orthogonality: ε(λ)·ε*(λ') = 0 for λ ≠ λ'
        for lam1 in 0..3 {
            for lam2 in 0..3 {
                if lam1 == lam2 {
                    continue;
                }
                let ec = conjugate(&eps[lam2]);
                let dot = dot_cc(&eps[lam1], &ec);
                assert!(
                    dot.norm() < 1e-10,
                    "ε(λ={})·ε*(λ'={}) = {:?}, expected 0",
                    lam1,
                    lam2,
                    dot
                );
            }
        }
    }

    #[test]
    fn test_polarization_completeness() {
        // Completeness relation: Σ_λ ε^μ(λ) ε*^ν(λ) = -g^{μν} + k^μ k^ν / M²
        let mw: f64 = 80.377;
        let e: f64 = 100.0;
        let p = (e * e - mw * mw).sqrt();
        let cos_theta: f64 = 0.3;
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
        let k = FourMomentum::new(e, p * sin_theta, 0.0, p * cos_theta);
        let k_arr = k.as_array();

        let eps = polarization_vectors(&k, mw);
        let metric = [1.0, -1.0, -1.0, -1.0];

        for mu in 0..4 {
            for nu in 0..4 {
                let mut sum = Complex64::new(0.0, 0.0);
                for lam in 0..3 {
                    let ec = conjugate(&eps[lam]);
                    sum += eps[lam][mu] * ec[nu];
                }
                // -g^{μν} + k^μ k^ν / M²
                // g^{μν} = diag(1,-1,-1,-1), so -g^{μν} = diag(-1,1,1,1)
                let expected = if mu == nu { -metric[mu] } else { 0.0 }
                    + k_arr[mu] * k_arr[nu] / (mw * mw);
                assert!(
                    (sum - Complex64::new(expected, 0.0)).norm() < 1e-8,
                    "completeness fails at ({}, {}): got {:?}, expected {}",
                    mu,
                    nu,
                    sum,
                    expected
                );
            }
        }
    }
}
