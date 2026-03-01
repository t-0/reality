use super::FourMomentum;

/// Kallen triangle function: lambda(a, b, c) = a^2 + b^2 + c^2 - 2ab - 2ac - 2bc
pub fn kallen(a: f64, b: f64, c: f64) -> f64 {
    a * a + b * b + c * c - 2.0 * a * b - 2.0 * a * c - 2.0 * b * c
}

/// Mandelstam variables for 2 -> 2 scattering: p1 + p2 -> p3 + p4.
#[derive(Debug, Clone, Copy)]
pub struct Mandelstam {
    pub s: f64,
    pub t: f64,
    pub u: f64,
}

impl Mandelstam {
    /// Compute from four-momenta p1 + p2 -> p3 + p4.
    pub fn from_momenta(
        p1: &FourMomentum,
        p2: &FourMomentum,
        p3: &FourMomentum,
        p4: &FourMomentum,
    ) -> Self {
        let s = (*p1 + *p2).invariant_mass_sq();
        let t = (*p1 - *p3).invariant_mass_sq();
        let u = (*p1 - *p4).invariant_mass_sq();
        Self { s, t, u }
    }

    /// For massless 2->2: compute from s and cos(theta) in the CM frame.
    pub fn from_s_costheta_massless(s: f64, cos_theta: f64) -> Self {
        let t = -s / 2.0 * (1.0 - cos_theta);
        let u = -s / 2.0 * (1.0 + cos_theta);
        Self { s, t, u }
    }

    /// For massive 2->2: compute from s, cos(theta), and initial/final masses in CM frame.
    /// Assumes symmetric initial state masses (m_i) and allows different final state masses.
    pub fn from_s_costheta(
        s: f64,
        cos_theta: f64,
        m_initial: [f64; 2],
        m_final: [f64; 2],
    ) -> Self {
        let sqrt_s = s.sqrt();
        let m1 = m_initial[0];
        let m2 = m_initial[1];
        let m3 = m_final[0];
        let m4 = m_final[1];

        // CM momentum magnitudes
        let p_i = kallen(s, m1 * m1, m2 * m2).max(0.0).sqrt() / (2.0 * sqrt_s);
        let p_f = kallen(s, m3 * m3, m4 * m4).max(0.0).sqrt() / (2.0 * sqrt_s);

        // Energies in CM
        let e1 = (s + m1 * m1 - m2 * m2) / (2.0 * sqrt_s);
        let e3 = (s + m3 * m3 - m4 * m4) / (2.0 * sqrt_s);

        let t = m1 * m1 + m3 * m3 - 2.0 * e1 * e3 + 2.0 * p_i * p_f * cos_theta;
        let u = m1 * m1 + m4 * m4 - 2.0 * e1 * (sqrt_s - e3) - 2.0 * p_i * p_f * cos_theta;

        Self { s, t, u }
    }

    /// Check the Mandelstam relation: s + t + u = sum of masses squared.
    pub fn sum(&self) -> f64 {
        self.s + self.t + self.u
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kallen() {
        // kallen(s, 0, 0) = s^2
        assert!((kallen(100.0, 0.0, 0.0) - 10000.0).abs() < 1e-10);
        // kallen(4, 1, 1) = 16 + 1 + 1 - 8 - 8 - 2 = 0
        assert!(kallen(4.0, 1.0, 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_mandelstam_massless() {
        let s = 100.0;
        let cos_theta = 0.5;
        let m = Mandelstam::from_s_costheta_massless(s, cos_theta);
        // s + t + u = 0 for massless
        assert!((m.sum()).abs() < 1e-10);
        // t = -s/2 (1-cos) = -25, u = -s/2 (1+cos) = -75
        assert!((m.t - (-25.0)).abs() < 1e-10);
        assert!((m.u - (-75.0)).abs() < 1e-10);
    }

    #[test]
    fn test_mandelstam_from_momenta() {
        // e+e- -> mu+mu- at sqrt(s) = 10 GeV, all massless
        let p1 = FourMomentum::new(5.0, 0.0, 0.0, 5.0);
        let p2 = FourMomentum::new(5.0, 0.0, 0.0, -5.0);
        let p3 = FourMomentum::new(5.0, 5.0, 0.0, 0.0);
        let p4 = FourMomentum::new(5.0, -5.0, 0.0, 0.0);
        let m = Mandelstam::from_momenta(&p1, &p2, &p3, &p4);
        assert!((m.s - 100.0).abs() < 1e-10);
        assert!((m.sum()).abs() < 1e-10);
    }
}
