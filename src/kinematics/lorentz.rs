use super::FourMomentum;

/// Boost a four-vector to the rest frame of the given reference momentum.
pub fn boost_to_rest(p: &FourMomentum, reference: &FourMomentum) -> FourMomentum {
    let beta = reference.beta();
    p.boost(&beta)
}

/// Boost a four-vector from the rest frame of the reference momentum back to the lab frame.
pub fn boost_from_rest(p: &FourMomentum, reference: &FourMomentum) -> FourMomentum {
    let beta = reference.beta();
    let neg_beta = [-beta[0], -beta[1], -beta[2]];
    p.boost(&neg_beta)
}

/// Rotate a four-momentum by angle phi around the z-axis.
pub fn rotate_z(p: &FourMomentum, phi: f64) -> FourMomentum {
    let cos_phi = phi.cos();
    let sin_phi = phi.sin();
    FourMomentum::new(
        p.e,
        p.px * cos_phi - p.py * sin_phi,
        p.px * sin_phi + p.py * cos_phi,
        p.pz,
    )
}

/// Rotate a four-momentum by angle theta around the y-axis.
pub fn rotate_y(p: &FourMomentum, theta: f64) -> FourMomentum {
    let cos_t = theta.cos();
    let sin_t = theta.sin();
    FourMomentum::new(
        p.e,
        p.px * cos_t + p.pz * sin_t,
        p.py,
        -p.px * sin_t + p.pz * cos_t,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_boost_to_rest_and_back() {
        let total = FourMomentum::new(15.0, 1.0, 2.0, 3.0);
        let p = FourMomentum::new(5.0, 0.5, 1.0, 1.5);
        let rest = boost_to_rest(&p, &total);
        let back = boost_from_rest(&rest, &total);
        assert!((back.e - p.e).abs() < 1e-9);
        assert!((back.px - p.px).abs() < 1e-9);
        assert!((back.py - p.py).abs() < 1e-9);
        assert!((back.pz - p.pz).abs() < 1e-9);
    }

    #[test]
    fn test_rotate_z_full_circle() {
        let p = FourMomentum::new(5.0, 1.0, 0.0, 3.0);
        let rotated = rotate_z(&p, 2.0 * std::f64::consts::PI);
        assert!((rotated.px - p.px).abs() < 1e-10);
        assert!((rotated.py - p.py).abs() < 1e-10);
    }

    #[test]
    fn test_rotate_y_preserves_energy() {
        let p = FourMomentum::new(5.0, 1.0, 2.0, 3.0);
        let rotated = rotate_y(&p, 0.7);
        assert!((rotated.e - p.e).abs() < 1e-10);
        let p_mag = p.three_momentum_mag();
        let r_mag = rotated.three_momentum_mag();
        assert!((p_mag - r_mag).abs() < 1e-10);
    }
}
