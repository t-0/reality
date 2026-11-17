use std::f64::consts::PI;

use super::field::Field;
use crate::couplings::constants::M_Z;
use crate::couplings::unification::Q_MAX;

/// Gauge boson self-contribution -11/3 C2(G) for U(1), SU(2), SU(3).
pub const GAUGE_BETA: [f64; 3] = [0.0, -22.0 / 3.0, -11.0];

/// One-loop beta coefficients (b1, b2, b3) from a list of active fields:
/// b = -11/3 C2(G) + 2/3 sum_Weyl T(R) + 1/3 sum_complex_scalar T(R).
pub fn beta_coefficients<'a>(fields: impl IntoIterator<Item = &'a Field>) -> [f64; 3] {
    let mut b = GAUGE_BETA;
    for f in fields {
        let c = f.beta_contribution();
        for i in 0..3 {
            b[i] += c[i];
        }
    }
    b
}

/// A range of scales over which the beta coefficients are constant.
#[derive(Debug, Clone)]
pub struct Segment {
    pub lo: f64,
    pub hi: f64,
    pub b: [f64; 3],
}

/// One-loop running with sharp mass thresholds between M_Z and Q_MAX.
#[derive(Debug, Clone)]
pub struct Running {
    pub segments: Vec<Segment>,
}

impl Running {
    pub fn new(fields: &[Field]) -> Self {
        let mut edges: Vec<f64> = fields
            .iter()
            .map(|f| f.mass)
            .filter(|&m| m > M_Z && m < Q_MAX)
            .collect();
        edges.push(M_Z);
        edges.push(Q_MAX);
        edges.sort_by(|a, b| a.partial_cmp(b).unwrap());
        edges.dedup_by(|a, b| (*a / *b - 1.0).abs() < 1e-9);

        let segments = edges
            .windows(2)
            .map(|w| Segment {
                lo: w[0],
                hi: w[1],
                b: beta_coefficients(fields.iter().filter(|f| f.mass <= w[0] * (1.0 + 1e-9))),
            })
            .collect();
        Self { segments }
    }

    /// Accumulated running Delta_i(q) = sum b_i/(2 pi) ln(...) from M_Z to q,
    /// so that 1/alpha_i(q) = 1/alpha_i(M_Z) - Delta_i(q).
    pub fn delta(&self, q: f64) -> [f64; 3] {
        let mut d = [0.0; 3];
        for s in &self.segments {
            if q <= s.lo {
                break;
            }
            let l = (q.min(s.hi) / s.lo).ln() / (2.0 * PI);
            for i in 0..3 {
                d[i] += s.b[i] * l;
            }
        }
        d
    }

    pub fn inverse_at(&self, inv_mz: [f64; 3], q: f64) -> [f64; 3] {
        let d = self.delta(q);
        [inv_mz[0] - d[0], inv_mz[1] - d[1], inv_mz[2] - d[2]]
    }

    /// Lowest scale between M_Z and Q_MAX where couplings i and j meet.
    /// `shift` is added to each 1/alpha before comparing (GUT threshold matching).
    /// At one loop the difference is linear in ln Q within each segment, so the
    /// crossing is solved exactly segment by segment.
    pub fn crossing(&self, inv_mz: [f64; 3], i: usize, j: usize, shift: [f64; 3]) -> Option<f64> {
        let mut d = (inv_mz[i] + shift[i]) - (inv_mz[j] + shift[j]);
        for s in &self.segments {
            if d == 0.0 {
                return Some(s.lo);
            }
            let slope = -(s.b[i] - s.b[j]) / (2.0 * PI);
            let d_hi = d + slope * (s.hi / s.lo).ln();
            if d_hi == 0.0 || d.signum() != d_hi.signum() {
                return Some(s.lo * (-d / slope).exp());
            }
            d = d_hi;
        }
        None
    }

    /// Lowest scale where any coupling becomes non-perturbative (1/alpha <= 0).
    pub fn landau_pole(&self, inv_mz: [f64; 3]) -> Option<f64> {
        for s in &self.segments {
            let a_hi = self.inverse_at(inv_mz, s.hi);
            if a_hi.iter().any(|&x| x <= 0.0) {
                // 1/alpha is linear in ln Q within a segment.
                let a_lo = self.inverse_at(inv_mz, s.lo);
                let mut first = s.hi;
                for i in 0..3 {
                    if a_hi[i] <= 0.0 {
                        let t = a_lo[i] / (a_lo[i] - a_hi[i]);
                        first = first.min(s.lo * (s.hi / s.lo).powf(t));
                    }
                }
                return Some(first);
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::couplings::unification::{inverse_couplings, inverse_couplings_mz, B_MSSM, B_SM};
    use crate::model::field::{Kind, Su3};
    use crate::model::sm::standard_model;

    /// MSSM superpartners (gauginos, higgsinos, second Higgs, sfermions).
    fn superpartners(mass: f64) -> Vec<Field> {
        vec![
            Field::new("gluino", Kind::Weyl, Su3::Octet, 1, 0.0, 1),
            Field::new("wino", Kind::Weyl, Su3::Singlet, 3, 0.0, 1),
            Field::new("bino", Kind::Weyl, Su3::Singlet, 1, 0.0, 1),
            Field::new("higgsinos", Kind::Dirac, Su3::Singlet, 2, 0.5, 1),
            Field::new("H_d", Kind::ComplexScalar, Su3::Singlet, 2, -0.5, 1),
            Field::new("squark_Q", Kind::ComplexScalar, Su3::Triplet, 2, 1.0 / 6.0, 3),
            Field::new("squark_u", Kind::ComplexScalar, Su3::AntiTriplet, 1, -2.0 / 3.0, 3),
            Field::new("squark_d", Kind::ComplexScalar, Su3::AntiTriplet, 1, 1.0 / 3.0, 3),
            Field::new("slepton_L", Kind::ComplexScalar, Su3::Singlet, 2, -0.5, 3),
            Field::new("slepton_e", Kind::ComplexScalar, Su3::Singlet, 1, 1.0, 3),
        ]
        .into_iter()
        .map(|f| f.with_mass(mass))
        .collect()
    }

    fn close(a: [f64; 3], b: [f64; 3]) -> bool {
        (0..3).all(|i| (a[i] - b[i]).abs() < 1e-12)
    }

    #[test]
    fn test_sm_beta_derived_from_content() {
        assert!(close(beta_coefficients(&standard_model()), B_SM));
    }

    #[test]
    fn test_mssm_beta_derived_from_content() {
        let mut fields = standard_model();
        fields.extend(superpartners(M_Z));
        assert!(close(beta_coefficients(&fields), B_MSSM));
    }

    #[test]
    fn test_threshold_matches_reference() {
        let mut fields = standard_model();
        fields.extend(superpartners(1000.0));
        let running = Running::new(&fields);
        assert_eq!(running.segments.len(), 2);
        let inv = inverse_couplings_mz();
        for q in [500.0, 1e6, 1e16] {
            let a = running.inverse_at(inv, q);
            let b = inverse_couplings(q, Some(1000.0));
            assert!(close_tol(a, b, 1e-9), "q = {}: {:?} vs {:?}", q, a, b);
        }
    }

    fn close_tol(a: [f64; 3], b: [f64; 3], tol: f64) -> bool {
        (0..3).all(|i| (a[i] - b[i]).abs() < tol)
    }

    #[test]
    fn test_landau_pole() {
        let inv = inverse_couplings_mz();
        assert!(Running::new(&standard_model()).landau_pole(inv).is_none());
        // Many charged scalars drive alpha_1 to a pole well below the Planck scale.
        let mut fields = standard_model();
        fields.push(Field::new("junk", Kind::ComplexScalar, Su3::Singlet, 1, 2.0, 40));
        let pole = Running::new(&fields).landau_pole(inv).unwrap();
        assert!(pole > M_Z && pole < 1e10, "pole at {:e}", pole);
        let a = Running::new(&fields).inverse_at(inv, pole);
        assert!(a[0].abs() < 1e-6);
    }
}
