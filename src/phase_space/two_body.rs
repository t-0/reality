use std::f64::consts::PI;

use rand::Rng;

use crate::kinematics::FourMomentum;

/// Kallen triangle function: lambda(a, b, c) = a^2 + b^2 + c^2 - 2ab - 2ac - 2bc
pub fn kallen_lambda(a: f64, b: f64, c: f64) -> f64 {
    a * a + b * b + c * c - 2.0 * a * b - 2.0 * a * c - 2.0 * b * c
}

/// Center-of-mass momentum for a two-body decay/scattering.
/// |p| = sqrt(lambda(s, m1^2, m2^2)) / (2 * sqrt(s))
pub fn two_body_momentum(sqrt_s: f64, m1: f64, m2: f64) -> f64 {
    let s = sqrt_s * sqrt_s;
    let lam = kallen_lambda(s, m1 * m1, m2 * m2);
    if lam <= 0.0 {
        return 0.0;
    }
    lam.sqrt() / (2.0 * sqrt_s)
}

/// Two-body phase space volume: PS_2 = |p_f| / (8 * pi * sqrt(s))
pub fn two_body_phase_space_volume(sqrt_s: f64, m1: f64, m2: f64) -> f64 {
    let p = two_body_momentum(sqrt_s, m1, m2);
    p / (8.0 * PI * sqrt_s)
}

/// Generate a two-body final state isotropically in the CM frame.
/// Returns (p3, p4) four-momenta for particles with masses m1, m2.
pub fn generate_two_body(
    rng: &mut impl Rng,
    sqrt_s: f64,
    m1: f64,
    m2: f64,
) -> Option<(FourMomentum, FourMomentum)> {
    let p = two_body_momentum(sqrt_s, m1, m2);
    if p <= 0.0 {
        return None;
    }

    let s = sqrt_s * sqrt_s;
    let e1 = (s + m1 * m1 - m2 * m2) / (2.0 * sqrt_s);
    let e2 = (s + m2 * m2 - m1 * m1) / (2.0 * sqrt_s);

    // Isotropic direction
    let cos_theta: f64 = rng.random_range(-1.0..1.0);
    let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
    let phi: f64 = rng.random_range(0.0..(2.0 * PI));

    let px = p * sin_theta * phi.cos();
    let py = p * sin_theta * phi.sin();
    let pz = p * cos_theta;

    let p3 = FourMomentum::new(e1, px, py, pz);
    let p4 = FourMomentum::new(e2, -px, -py, -pz);

    Some((p3, p4))
}

/// Generate a two-body final state at a specific cos(theta) and phi.
pub fn two_body_at_angle(
    sqrt_s: f64,
    m1: f64,
    m2: f64,
    cos_theta: f64,
    phi: f64,
) -> Option<(FourMomentum, FourMomentum)> {
    let p = two_body_momentum(sqrt_s, m1, m2);
    if p <= 0.0 {
        return None;
    }

    let s = sqrt_s * sqrt_s;
    let e1 = (s + m1 * m1 - m2 * m2) / (2.0 * sqrt_s);
    let e2 = (s + m2 * m2 - m1 * m1) / (2.0 * sqrt_s);

    let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
    let px = p * sin_theta * phi.cos();
    let py = p * sin_theta * phi.sin();
    let pz = p * cos_theta;

    let p3 = FourMomentum::new(e1, px, py, pz);
    let p4 = FourMomentum::new(e2, -px, -py, -pz);

    Some((p3, p4))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_two_body_momentum_massless() {
        let p = two_body_momentum(10.0, 0.0, 0.0);
        assert!((p - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_two_body_momentum_threshold() {
        // At threshold: sqrt(s) = m1 + m2, |p| = 0
        let p = two_body_momentum(2.0, 1.0, 1.0);
        assert!(p.abs() < 1e-10);
    }

    #[test]
    fn test_two_body_momentum_below_threshold() {
        let p = two_body_momentum(1.0, 1.0, 1.0);
        assert_eq!(p, 0.0);
    }

    #[test]
    fn test_generate_conserves_momentum() {
        let mut rng = rand::rng();
        let sqrt_s = 200.0;
        let m1 = 80.377;
        let m2 = 80.377;

        let (p3, p4) = generate_two_body(&mut rng, sqrt_s, m1, m2).unwrap();
        let total = p3 + p4;

        assert!((total.e - sqrt_s).abs() < 1e-10);
        assert!(total.px.abs() < 1e-10);
        assert!(total.py.abs() < 1e-10);
        assert!(total.pz.abs() < 1e-10);
        assert!((p3.invariant_mass() - m1).abs() < 1e-6);
        assert!((p4.invariant_mass() - m2).abs() < 1e-6);
    }
}
