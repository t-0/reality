use num_complex::Complex64;
use std::ops::{Add, Mul, Sub};

use crate::kinematics::FourMomentum;

const ZERO: Complex64 = Complex64::new(0.0, 0.0);
const ONE: Complex64 = Complex64::new(1.0, 0.0);
const NEG_ONE: Complex64 = Complex64::new(-1.0, 0.0);
const I: Complex64 = Complex64::new(0.0, 1.0);
const NEG_I: Complex64 = Complex64::new(0.0, -1.0);

/// A 4x4 complex matrix used for Dirac algebra.
/// Stored as row-major: data[row][col].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DiracMatrix {
    pub data: [[Complex64; 4]; 4],
}

impl DiracMatrix {
    pub const ZERO: Self = Self {
        data: [[ZERO; 4]; 4],
    };

    pub const IDENTITY: Self = Self {
        data: [
            [ONE, ZERO, ZERO, ZERO],
            [ZERO, ONE, ZERO, ZERO],
            [ZERO, ZERO, ONE, ZERO],
            [ZERO, ZERO, ZERO, ONE],
        ],
    };

    /// Dirac representation gamma matrices.
    /// gamma^0 = diag(1,1,-1,-1)
    pub fn gamma0() -> Self {
        Self {
            data: [
                [ONE, ZERO, ZERO, ZERO],
                [ZERO, ONE, ZERO, ZERO],
                [ZERO, ZERO, NEG_ONE, ZERO],
                [ZERO, ZERO, ZERO, NEG_ONE],
            ],
        }
    }

    /// gamma^1
    pub fn gamma1() -> Self {
        Self {
            data: [
                [ZERO, ZERO, ZERO, ONE],
                [ZERO, ZERO, ONE, ZERO],
                [ZERO, NEG_ONE, ZERO, ZERO],
                [NEG_ONE, ZERO, ZERO, ZERO],
            ],
        }
    }

    /// gamma^2
    pub fn gamma2() -> Self {
        Self {
            data: [
                [ZERO, ZERO, ZERO, NEG_I],
                [ZERO, ZERO, I, ZERO],
                [ZERO, I, ZERO, ZERO],
                [NEG_I, ZERO, ZERO, ZERO],
            ],
        }
    }

    /// gamma^3
    pub fn gamma3() -> Self {
        Self {
            data: [
                [ZERO, ZERO, ONE, ZERO],
                [ZERO, ZERO, ZERO, NEG_ONE],
                [NEG_ONE, ZERO, ZERO, ZERO],
                [ZERO, ONE, ZERO, ZERO],
            ],
        }
    }

    /// Return gamma^mu for mu = 0, 1, 2, 3.
    pub fn gamma(mu: usize) -> Self {
        match mu {
            0 => Self::gamma0(),
            1 => Self::gamma1(),
            2 => Self::gamma2(),
            3 => Self::gamma3(),
            _ => panic!("gamma index must be 0..3, got {}", mu),
        }
    }

    /// gamma^5 = i * gamma^0 * gamma^1 * gamma^2 * gamma^3
    pub fn gamma5() -> Self {
        // In Dirac representation:
        // gamma5 = [[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]
        Self {
            data: [
                [ZERO, ZERO, ONE, ZERO],
                [ZERO, ZERO, ZERO, ONE],
                [ONE, ZERO, ZERO, ZERO],
                [ZERO, ONE, ZERO, ZERO],
            ],
        }
    }

    /// Left-handed projection operator: P_L = (1 - gamma5) / 2
    pub fn proj_left() -> Self {
        let half = Complex64::new(0.5, 0.0);
        (Self::IDENTITY - Self::gamma5()).scale(half)
    }

    /// Right-handed projection operator: P_R = (1 + gamma5) / 2
    pub fn proj_right() -> Self {
        let half = Complex64::new(0.5, 0.0);
        (Self::IDENTITY + Self::gamma5()).scale(half)
    }

    /// Feynman slash: p_slash = p_mu * gamma^mu = E*gamma0 - px*gamma1 - py*gamma2 - pz*gamma3
    pub fn slash(p: &FourMomentum) -> Self {
        let mut result = Self::ZERO;
        let g = [Self::gamma0(), Self::gamma1(), Self::gamma2(), Self::gamma3()];
        let metric = [1.0, -1.0, -1.0, -1.0];
        let components = [p.e, p.px, p.py, p.pz];

        for mu in 0..4 {
            // p_mu = g_{mu mu} * p^mu, and gamma^mu is upper index
            // p_slash = p_mu gamma^mu = p^0 gamma^0 - p^1 gamma^1 - p^2 gamma^2 - p^3 gamma^3
            // But we store p as contravariant (E, px, py, pz), so:
            // p_slash = E gamma^0 - px gamma^1 - py gamma^2 - pz gamma^3
            // Actually p_slash = p_mu gamma^mu where p_mu = (E, -px, -py, -pz)
            // So p_slash = E*gamma0 + (-px)*gamma1 + (-py)*gamma2 + (-pz)*gamma3
            // Which is the same since gamma_mu = g_{mu nu} gamma^nu, and the metric
            // absorbs the signs. Let's use p_slash = p^mu gamma_mu where gamma_mu raises:
            // p_slash = p^mu g_{mu nu} gamma^nu = E gamma^0 - px gamma^1 - py gamma^2 - pz gamma^3
            let coeff = Complex64::new(metric[mu] * components[mu], 0.0);
            for i in 0..4 {
                for j in 0..4 {
                    result.data[i][j] += coeff * g[mu].data[i][j];
                }
            }
        }
        result
    }

    /// Feynman slash for a complex 4-vector: ε̸ = ε_μ γ^μ = ε^0 γ^0 - ε^1 γ^1 - ε^2 γ^2 - ε^3 γ^3
    /// Input components are contravariant; the metric is applied internally.
    pub fn slash_complex(eps: &[Complex64; 4]) -> Self {
        let mut result = Self::ZERO;
        let g = [Self::gamma0(), Self::gamma1(), Self::gamma2(), Self::gamma3()];
        let metric = [1.0, -1.0, -1.0, -1.0];

        for mu in 0..4 {
            let coeff = Complex64::new(metric[mu], 0.0) * eps[mu];
            for i in 0..4 {
                for j in 0..4 {
                    result.data[i][j] += coeff * g[mu].data[i][j];
                }
            }
        }
        result
    }

    /// Trace of this matrix: Tr[M] = sum_i M[i][i]
    pub fn trace(&self) -> Complex64 {
        self.data[0][0] + self.data[1][1] + self.data[2][2] + self.data[3][3]
    }

    /// Matrix multiplication: self * other
    pub fn multiply(&self, other: &Self) -> Self {
        let mut result = Self::ZERO;
        for i in 0..4 {
            for j in 0..4 {
                let mut sum = ZERO;
                for k in 0..4 {
                    sum += self.data[i][k] * other.data[k][j];
                }
                result.data[i][j] = sum;
            }
        }
        result
    }

    /// Hermitian conjugate (dagger): M^dagger[i][j] = conj(M[j][i])
    pub fn dagger(&self) -> Self {
        let mut result = Self::ZERO;
        for i in 0..4 {
            for j in 0..4 {
                result.data[i][j] = self.data[j][i].conj();
            }
        }
        result
    }

    /// Dirac adjoint: M_bar = gamma0 * M^dagger * gamma0
    pub fn adjoint(&self) -> Self {
        Self::gamma0().multiply(&self.dagger()).multiply(&Self::gamma0())
    }

    /// Scale all elements by a complex constant.
    pub fn scale(&self, c: Complex64) -> Self {
        let mut result = *self;
        for i in 0..4 {
            for j in 0..4 {
                result.data[i][j] *= c;
            }
        }
        result
    }

    /// Scale by a real constant.
    pub fn scale_real(&self, c: f64) -> Self {
        self.scale(Complex64::new(c, 0.0))
    }
}

impl Add for DiracMatrix {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let mut result = Self::ZERO;
        for i in 0..4 {
            for j in 0..4 {
                result.data[i][j] = self.data[i][j] + rhs.data[i][j];
            }
        }
        result
    }
}

impl Sub for DiracMatrix {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let mut result = Self::ZERO;
        for i in 0..4 {
            for j in 0..4 {
                result.data[i][j] = self.data[i][j] - rhs.data[i][j];
            }
        }
        result
    }
}

impl Mul for DiracMatrix {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        self.multiply(&rhs)
    }
}

impl Mul<Complex64> for DiracMatrix {
    type Output = Self;
    fn mul(self, rhs: Complex64) -> Self {
        self.scale(rhs)
    }
}

impl Mul<f64> for DiracMatrix {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        self.scale_real(rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_clifford_algebra() {
        // {gamma^mu, gamma^nu} = 2 * g^{mu nu} * I
        let metric = [1.0, -1.0, -1.0, -1.0];
        for mu in 0..4 {
            for nu in 0..4 {
                let gm = DiracMatrix::gamma(mu);
                let gn = DiracMatrix::gamma(nu);
                let anticomm = gm * gn + gn * gm;
                let expected = if mu == nu {
                    DiracMatrix::IDENTITY.scale_real(2.0 * metric[mu])
                } else {
                    DiracMatrix::ZERO
                };
                for i in 0..4 {
                    for j in 0..4 {
                        assert!(
                            (anticomm.data[i][j] - expected.data[i][j]).norm() < 1e-10,
                            "Clifford failed for mu={}, nu={}, i={}, j={}",
                            mu, nu, i, j
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn test_gamma5_squared() {
        let g5 = DiracMatrix::gamma5();
        let g5sq = g5 * g5;
        for i in 0..4 {
            for j in 0..4 {
                let expected = if i == j { ONE } else { ZERO };
                assert!(
                    (g5sq.data[i][j] - expected).norm() < 1e-10,
                    "gamma5^2 != I at ({}, {})",
                    i, j
                );
            }
        }
    }

    #[test]
    fn test_gamma5_anticommutes() {
        // {gamma5, gamma^mu} = 0 for mu = 0..3
        let g5 = DiracMatrix::gamma5();
        for mu in 0..4 {
            let gm = DiracMatrix::gamma(mu);
            let anticomm = g5 * gm + gm * g5;
            for i in 0..4 {
                for j in 0..4 {
                    assert!(
                        anticomm.data[i][j].norm() < 1e-10,
                        "{{gamma5, gamma{}}} != 0 at ({}, {})",
                        mu, i, j
                    );
                }
            }
        }
    }

    #[test]
    fn test_trace_identity() {
        assert!((DiracMatrix::IDENTITY.trace() - Complex64::new(4.0, 0.0)).norm() < 1e-10);
    }

    #[test]
    fn test_trace_gamma() {
        // Tr[gamma^mu] = 0
        for mu in 0..4 {
            assert!(DiracMatrix::gamma(mu).trace().norm() < 1e-10);
        }
    }

    #[test]
    fn test_trace_gamma5() {
        assert!(DiracMatrix::gamma5().trace().norm() < 1e-10);
    }

    #[test]
    fn test_trace_two_gammas() {
        // Tr[gamma^mu gamma^nu] = 4 * g^{mu nu}
        let metric = [1.0, -1.0, -1.0, -1.0];
        for mu in 0..4 {
            for nu in 0..4 {
                let tr = (DiracMatrix::gamma(mu) * DiracMatrix::gamma(nu)).trace();
                let expected = if mu == nu {
                    Complex64::new(4.0 * metric[mu], 0.0)
                } else {
                    ZERO
                };
                assert!(
                    (tr - expected).norm() < 1e-10,
                    "Tr[gamma{} gamma{}] = {:?}, expected {:?}",
                    mu, nu, tr, expected
                );
            }
        }
    }

    #[test]
    fn test_slash_massless() {
        // For a massless particle along z: p = (E, 0, 0, E)
        // p_slash^2 = p^2 * I = 0
        let p = FourMomentum::new(5.0, 0.0, 0.0, 5.0);
        let ps = DiracMatrix::slash(&p);
        let ps2 = ps * ps;
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    ps2.data[i][j].norm() < 1e-10,
                    "slash(p)^2 != 0 for massless at ({}, {})",
                    i, j
                );
            }
        }
    }

    #[test]
    fn test_slash_massive() {
        // For a massive particle: p_slash^2 = m^2 * I
        let m = 1.0;
        let e = (m * m + 9.0_f64).sqrt();
        let p = FourMomentum::new(e, 0.0, 0.0, 3.0);
        let ps = DiracMatrix::slash(&p);
        let ps2 = ps * ps;
        for i in 0..4 {
            for j in 0..4 {
                let expected = if i == j {
                    Complex64::new(m * m, 0.0)
                } else {
                    ZERO
                };
                assert!(
                    (ps2.data[i][j] - expected).norm() < 1e-10,
                    "slash(p)^2 != m^2 * I at ({}, {})",
                    i, j
                );
            }
        }
    }

    #[test]
    fn test_projectors() {
        let pl = DiracMatrix::proj_left();
        let pr = DiracMatrix::proj_right();

        // PL + PR = I
        let sum = pl + pr;
        for i in 0..4 {
            for j in 0..4 {
                let expected = if i == j { ONE } else { ZERO };
                assert!((sum.data[i][j] - expected).norm() < 1e-10);
            }
        }

        // PL * PL = PL
        let pl2 = pl * pl;
        for i in 0..4 {
            for j in 0..4 {
                assert!((pl2.data[i][j] - pl.data[i][j]).norm() < 1e-10);
            }
        }

        // PL * PR = 0
        let plpr = pl * pr;
        for i in 0..4 {
            for j in 0..4 {
                assert!(plpr.data[i][j].norm() < 1e-10);
            }
        }
    }
}
