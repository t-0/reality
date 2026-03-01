use std::ops::{Add, Mul, Neg, Sub};

/// A Lorentz four-vector in the (E, px, py, pz) convention.
/// Metric signature: (+, -, -, -).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FourMomentum {
    pub e: f64,
    pub px: f64,
    pub py: f64,
    pub pz: f64,
}

impl FourMomentum {
    pub fn new(e: f64, px: f64, py: f64, pz: f64) -> Self {
        Self { e, px, py, pz }
    }

    /// A particle at rest with given mass.
    pub fn at_rest(mass: f64) -> Self {
        Self::new(mass, 0.0, 0.0, 0.0)
    }

    /// Zero four-vector.
    pub fn zero() -> Self {
        Self::new(0.0, 0.0, 0.0, 0.0)
    }

    /// Minkowski inner product: p . q = p0*q0 - p1*q1 - p2*q2 - p3*q3
    pub fn dot(&self, other: &Self) -> f64 {
        self.e * other.e - self.px * other.px - self.py * other.py - self.pz * other.pz
    }

    /// Invariant mass squared: p^2 = E^2 - |p|^2
    pub fn invariant_mass_sq(&self) -> f64 {
        self.dot(self)
    }

    /// Invariant mass (returns 0 if mass_sq is slightly negative due to numerics).
    pub fn invariant_mass(&self) -> f64 {
        self.invariant_mass_sq().max(0.0).sqrt()
    }

    /// Magnitude of the three-momentum.
    pub fn three_momentum_mag(&self) -> f64 {
        (self.px * self.px + self.py * self.py + self.pz * self.pz).sqrt()
    }

    /// Transverse momentum pT = sqrt(px^2 + py^2).
    pub fn transverse_momentum(&self) -> f64 {
        (self.px * self.px + self.py * self.py).sqrt()
    }

    /// Rapidity: y = 0.5 * ln((E + pz) / (E - pz))
    pub fn rapidity(&self) -> f64 {
        0.5 * ((self.e + self.pz) / (self.e - self.pz)).ln()
    }

    /// Pseudorapidity: eta = -ln(tan(theta/2))
    pub fn pseudorapidity(&self) -> f64 {
        let p = self.three_momentum_mag();
        0.5 * ((p + self.pz) / (p - self.pz)).ln()
    }

    /// Spatial components as an array.
    pub fn spatial(&self) -> [f64; 3] {
        [self.px, self.py, self.pz]
    }

    /// Velocity vector beta = p/E.
    pub fn beta(&self) -> [f64; 3] {
        [self.px / self.e, self.py / self.e, self.pz / self.e]
    }

    /// Lorentz boost by velocity beta = (bx, by, bz).
    /// Transforms this four-vector to a frame moving with velocity -beta.
    pub fn boost(&self, beta: &[f64; 3]) -> Self {
        let bx = beta[0];
        let by = beta[1];
        let bz = beta[2];
        let b2 = bx * bx + by * by + bz * bz;

        if b2 < 1e-30 {
            return *self;
        }

        let gamma = 1.0 / (1.0 - b2).sqrt();
        let bp = bx * self.px + by * self.py + bz * self.pz;
        let gamma2 = if b2 > 0.0 {
            (gamma - 1.0) / b2
        } else {
            0.0
        };

        Self {
            e: gamma * (self.e - bp),
            px: self.px + (gamma2 * bp - gamma * self.e) * bx,
            py: self.py + (gamma2 * bp - gamma * self.e) * by,
            pz: self.pz + (gamma2 * bp - gamma * self.e) * bz,
        }
    }

    /// As a 4-element array [E, px, py, pz].
    pub fn as_array(&self) -> [f64; 4] {
        [self.e, self.px, self.py, self.pz]
    }
}

impl Add for FourMomentum {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            e: self.e + rhs.e,
            px: self.px + rhs.px,
            py: self.py + rhs.py,
            pz: self.pz + rhs.pz,
        }
    }
}

impl Sub for FourMomentum {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            e: self.e - rhs.e,
            px: self.px - rhs.px,
            py: self.py - rhs.py,
            pz: self.pz - rhs.pz,
        }
    }
}

impl Neg for FourMomentum {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            e: -self.e,
            px: -self.px,
            py: -self.py,
            pz: -self.pz,
        }
    }
}

impl Mul<f64> for FourMomentum {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self {
            e: self.e * rhs,
            px: self.px * rhs,
            py: self.py * rhs,
            pz: self.pz * rhs,
        }
    }
}

impl Mul<FourMomentum> for f64 {
    type Output = FourMomentum;
    fn mul(self, rhs: FourMomentum) -> FourMomentum {
        rhs * self
    }
}

impl Add for &FourMomentum {
    type Output = FourMomentum;
    fn add(self, rhs: Self) -> FourMomentum {
        *self + *rhs
    }
}

impl Sub for &FourMomentum {
    type Output = FourMomentum;
    fn sub(self, rhs: Self) -> FourMomentum {
        *self - *rhs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_invariant_mass_at_rest() {
        let p = FourMomentum::at_rest(0.511e-3); // electron at rest
        assert!((p.invariant_mass() - 0.511e-3).abs() < 1e-10);
    }

    #[test]
    fn test_invariant_mass_moving() {
        // Particle with mass 1 GeV, pz = 3 GeV => E = sqrt(1+9) = sqrt(10)
        let e = (1.0_f64 + 9.0).sqrt();
        let p = FourMomentum::new(e, 0.0, 0.0, 3.0);
        assert!((p.invariant_mass() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_minkowski_dot() {
        let p1 = FourMomentum::new(5.0, 1.0, 2.0, 3.0);
        let p2 = FourMomentum::new(4.0, 0.5, 1.0, 1.5);
        // dot = 5*4 - 1*0.5 - 2*1 - 3*1.5 = 20 - 0.5 - 2 - 4.5 = 13
        assert!((p1.dot(&p2) - 13.0).abs() < 1e-10);
    }

    #[test]
    fn test_massless_photon() {
        let p = FourMomentum::new(10.0, 0.0, 0.0, 10.0);
        assert!(p.invariant_mass_sq().abs() < 1e-10);
    }

    #[test]
    fn test_boost_round_trip() {
        let p = FourMomentum::new(10.0, 1.0, 2.0, 3.0);
        let beta = [0.3, 0.1, 0.2];
        let boosted = p.boost(&beta);
        let neg_beta = [-beta[0], -beta[1], -beta[2]];
        let recovered = boosted.boost(&neg_beta);
        assert!((recovered.e - p.e).abs() < 1e-10);
        assert!((recovered.px - p.px).abs() < 1e-10);
        assert!((recovered.py - p.py).abs() < 1e-10);
        assert!((recovered.pz - p.pz).abs() < 1e-10);
    }

    #[test]
    fn test_boost_to_rest() {
        // Particle with mass 1, moving along z with pz=3, E=sqrt(10)
        let e = 10.0_f64.sqrt();
        let p = FourMomentum::new(e, 0.0, 0.0, 3.0);
        let b = p.beta();
        let rest = p.boost(&b);
        assert!((rest.invariant_mass() - 1.0).abs() < 1e-10);
        assert!(rest.three_momentum_mag() < 1e-10);
    }

    #[test]
    fn test_addition() {
        let p1 = FourMomentum::new(5.0, 1.0, 0.0, 0.0);
        let p2 = FourMomentum::new(5.0, -1.0, 0.0, 0.0);
        let total = p1 + p2;
        assert!((total.e - 10.0).abs() < 1e-10);
        assert!(total.px.abs() < 1e-10);
        assert!((total.invariant_mass() - 10.0).abs() < 1e-10);
    }
}
