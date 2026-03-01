use num_complex::Complex64;

use super::gamma::DiracMatrix;
use crate::kinematics::FourMomentum;

const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const ONE: Complex64 = Complex64::new(1.0, 0.0);

/// Spin state for a fermion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SpinState {
    Up,
    Down,
}

/// A 4-component Dirac spinor (column vector).
#[derive(Debug, Clone, Copy)]
pub struct DiracSpinor {
    pub components: [Complex64; 4],
}

/// A Dirac adjoint spinor (row vector): ubar = u^dagger * gamma0.
#[derive(Debug, Clone, Copy)]
pub struct DiracSpinorBar {
    pub components: [Complex64; 4],
}

impl DiracSpinor {
    /// Positive-energy spinor u(p, s) in the Dirac representation.
    /// For a particle with four-momentum p and mass m.
    pub fn u(p: &FourMomentum, mass: f64, spin: SpinState) -> Self {
        let e = p.e;
        let px = p.px;
        let py = p.py;
        let pz = p.pz;

        let ep = (e + mass).max(0.0);
        let norm = (ep).sqrt();

        if norm < 1e-30 {
            return Self {
                components: [ZERO; 4],
            };
        }

        let chi = match spin {
            SpinState::Up => [ONE, ZERO],
            SpinState::Down => [ZERO, ONE],
        };

        // sigma . p / (E + m) applied to chi
        // sigma.p matrix = [[pz, px - i*py], [px + i*py, -pz]]
        let sigma_p_over_ep = if ep > 1e-30 {
            let inv = 1.0 / ep;
            let sp0 = Complex64::new(pz, 0.0) * chi[0]
                + Complex64::new(px, -py) * chi[1];
            let sp1 = Complex64::new(px, py) * chi[0]
                + Complex64::new(-pz, 0.0) * chi[1];
            [sp0 * inv, sp1 * inv]
        } else {
            [ZERO, ZERO]
        };

        Self {
            components: [
                norm * chi[0],
                norm * chi[1],
                norm * sigma_p_over_ep[0],
                norm * sigma_p_over_ep[1],
            ],
        }
    }

    /// Negative-energy spinor v(p, s) in the Dirac representation.
    /// For an antiparticle with four-momentum p and mass m.
    pub fn v(p: &FourMomentum, mass: f64, spin: SpinState) -> Self {
        let e = p.e;
        let px = p.px;
        let py = p.py;
        let pz = p.pz;

        let ep = (e + mass).max(0.0);
        let norm = (ep).sqrt();

        if norm < 1e-30 {
            return Self {
                components: [ZERO; 4],
            };
        }

        // For v-spinors, the spin labeling is reversed:
        // v(p, up) uses chi_down and v(p, down) uses chi_up
        let chi = match spin {
            SpinState::Up => [ZERO, ONE],
            SpinState::Down => [ONE, ZERO],
        };

        let sigma_p_over_ep = if ep > 1e-30 {
            let inv = 1.0 / ep;
            let sp0 =
                Complex64::new(pz, 0.0) * chi[0] + Complex64::new(px, -py) * chi[1];
            let sp1 =
                Complex64::new(px, py) * chi[0] + Complex64::new(-pz, 0.0) * chi[1];
            [sp0 * inv, sp1 * inv]
        } else {
            [ZERO, ZERO]
        };

        Self {
            components: [
                norm * sigma_p_over_ep[0],
                norm * sigma_p_over_ep[1],
                norm * chi[0],
                norm * chi[1],
            ],
        }
    }

    /// Dirac adjoint: ubar = u^dagger * gamma0
    pub fn bar(&self) -> DiracSpinorBar {
        let u_dag = [
            self.components[0].conj(),
            self.components[1].conj(),
            self.components[2].conj(),
            self.components[3].conj(),
        ];
        // gamma0 in Dirac rep: diag(1, 1, -1, -1)
        DiracSpinorBar {
            components: [u_dag[0], u_dag[1], -u_dag[2], -u_dag[3]],
        }
    }

    /// Spin sum for positive-energy spinors:
    /// sum_s u(p,s) ubar(p,s) = (p_slash + m * I)
    pub fn spin_sum_u(p: &FourMomentum, mass: f64) -> DiracMatrix {
        DiracMatrix::slash(p) + DiracMatrix::IDENTITY * mass
    }

    /// Spin sum for negative-energy spinors:
    /// sum_s v(p,s) vbar(p,s) = (p_slash - m * I)
    pub fn spin_sum_v(p: &FourMomentum, mass: f64) -> DiracMatrix {
        DiracMatrix::slash(p) - DiracMatrix::IDENTITY * mass
    }
}

impl DiracSpinorBar {
    /// Compute ubar . M . v (a scalar) where M is a 4x4 Dirac matrix.
    pub fn sandwich(&self, m: &DiracMatrix, v: &DiracSpinor) -> Complex64 {
        let mut result = ZERO;
        for i in 0..4 {
            for j in 0..4 {
                result += self.components[i] * m.data[i][j] * v.components[j];
            }
        }
        result
    }

    /// Compute the spinor product ubar . v (no matrix insertion).
    pub fn dot(&self, v: &DiracSpinor) -> Complex64 {
        let mut result = ZERO;
        for i in 0..4 {
            result += self.components[i] * v.components[i];
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spin_sum_u_at_rest() {
        let mass = 1.0;
        let p = FourMomentum::at_rest(mass);

        // Explicit spin sum
        let u_up = DiracSpinor::u(&p, mass, SpinState::Up);
        let u_down = DiracSpinor::u(&p, mass, SpinState::Down);

        let mut explicit = DiracMatrix::ZERO;
        for spin in [&u_up, &u_down] {
            let bar = spin.bar();
            for i in 0..4 {
                for j in 0..4 {
                    explicit.data[i][j] += spin.components[i] * bar.components[j];
                }
            }
        }

        let expected = DiracSpinor::spin_sum_u(&p, mass);
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (explicit.data[i][j] - expected.data[i][j]).norm() < 1e-10,
                    "spin sum u mismatch at ({}, {}): {:?} vs {:?}",
                    i, j, explicit.data[i][j], expected.data[i][j]
                );
            }
        }
    }

    #[test]
    fn test_spin_sum_u_moving() {
        // On-shell: E = sqrt(m^2 + |p|^2) with |p|^2 = 1+1+4 = 6
        let mass = 0.5;
        let e = (mass * mass + 6.0_f64).sqrt();
        let p = FourMomentum::new(e, 1.0, 1.0, 2.0);

        let u_up = DiracSpinor::u(&p, mass, SpinState::Up);
        let u_down = DiracSpinor::u(&p, mass, SpinState::Down);

        let mut explicit = DiracMatrix::ZERO;
        for spin in [&u_up, &u_down] {
            let bar = spin.bar();
            for i in 0..4 {
                for j in 0..4 {
                    explicit.data[i][j] += spin.components[i] * bar.components[j];
                }
            }
        }

        let expected = DiracSpinor::spin_sum_u(&p, mass);
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (explicit.data[i][j] - expected.data[i][j]).norm() < 1e-6,
                    "spin sum u moving mismatch at ({}, {})",
                    i, j
                );
            }
        }
    }

    #[test]
    fn test_ubar_u_normalization() {
        // ubar(p, s) . u(p, s) = 2m for each spin
        // On-shell: E = sqrt(m^2 + |p|^2) with |p|^2 = 1+1+4 = 6
        let mass = 1.5;
        let e = (mass * mass + 6.0_f64).sqrt();
        let p = FourMomentum::new(e, 1.0, 1.0, 2.0);

        for spin in [SpinState::Up, SpinState::Down] {
            let u = DiracSpinor::u(&p, mass, spin);
            let ubar = u.bar();
            let product = ubar.dot(&u);
            assert!(
                (product.re - 2.0 * mass).abs() < 1e-6,
                "ubar.u = {:?}, expected {}",
                product,
                2.0 * mass
            );
            assert!(product.im.abs() < 1e-6);
        }
    }
}
