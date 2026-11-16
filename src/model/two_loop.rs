//! Two-loop running of the gauge couplings with the top Yukawa, in t = ln Q:
//!
//!   d(1/alpha_i)/dt = -b_i/(2 pi) - sum_j b_ij alpha_j/(8 pi^2) + d_i y_t^2/(32 pi^3)
//!   d y_t/dt        = y_t/(16 pi^2) [a y_t^2 - sum_i c_i g_i^2]
//!
//! b_ij is derived from field content (Machacek & Vaughn, Nucl. Phys. B222 (1983) 83):
//!   b_ij = -34/3 C2(G_i)^2 delta_ij
//!        + sum_fermions  w [4 C2_j(R) + 20/3 C2(G_i) delta_ij] S2_i(R)
//!        + sum_scalars   w [4 C2_j(R) +  2/3 C2(G_i) delta_ij] S2_i(R)
//! with w = 1/2 per Weyl fermion and 1 per complex scalar.
//!
//! In supersymmetric ranges the gaugino-fermion-sfermion Yukawas (equal to the
//! gauge couplings) also contribute, and the superfield form is used instead:
//!   b_ij = -6 C2(G_i)^2 delta_ij + sum_chiral S2_i(Phi) [2 C2(G_i) delta_ij + 4 C2_j(Phi)]
//! A range counts as supersymmetric when it has one gaugino per gauge group and
//! every other Weyl fermion pairs with a complex scalar of the same charges.
//! Thresholds are sharp; threshold matching corrections are not included.

use std::f64::consts::PI;

use std::collections::HashMap;

use super::field::{Field, Kind, Su3};
use super::running::beta_coefficients;
use super::spec::TopYukawa;
use crate::couplings::constants::M_Z;
use crate::couplings::unification::Q_MAX;

/// Quadratic Casimir of the adjoint for U(1), SU(2), SU(3).
const GAUGE_CASIMIR: [f64; 3] = [0.0, 2.0, 3.0];

/// MS-bar top Yukawa at M_Z (from m_t(m_t) ~ 162.5 GeV run to M_Z).
pub const Y_TOP_MZ: f64 = 0.97;

/// Top Yukawa coefficients: (d_i in the gauge running, a, c_i in the Yukawa running).
#[derive(Debug, Clone, Copy)]
struct YukawaCoeffs {
    d: [f64; 3],
    a: f64,
    c: [f64; 3],
}

const SM_YUKAWA: YukawaCoeffs = YukawaCoeffs {
    d: [17.0 / 10.0, 3.0 / 2.0, 2.0],
    a: 9.0 / 2.0,
    c: [17.0 / 20.0, 9.0 / 4.0, 8.0],
};

const MSSM_YUKAWA: YukawaCoeffs = YukawaCoeffs {
    d: [26.0 / 5.0, 6.0, 4.0],
    a: 6.0,
    c: [13.0 / 15.0, 3.0, 16.0 / 3.0],
};

/// Two-loop gauge coefficient matrix b_ij from a list of active fields.
pub fn two_loop_matrix<'a>(fields: impl IntoIterator<Item = &'a Field>) -> [[f64; 3]; 3] {
    let mut b = [[0.0; 3]; 3];
    for i in 0..3 {
        b[i][i] = -34.0 / 3.0 * GAUGE_CASIMIR[i] * GAUGE_CASIMIR[i];
    }
    for f in fields {
        let w = f.kind.two_loop_weight();
        let gauge_term = if f.kind.is_fermion() { 20.0 / 3.0 } else { 2.0 / 3.0 };
        let s2 = f.dynkin();
        let c2 = f.casimirs();
        for i in 0..3 {
            for j in 0..3 {
                let diag = if i == j { gauge_term * GAUGE_CASIMIR[i] } else { 0.0 };
                b[i][j] += w * (4.0 * c2[j] + diag) * s2[i];
            }
        }
    }
    b
}

/// Gauge charges up to conjugation, as an exact key: (SU(3) code, SU(2) dim, 3600 Y).
fn charge_key(su3: Su3, su2: u32, y: f64) -> (u8, u32, i64) {
    let code = |r: Su3| match r {
        Su3::Singlet => 0,
        Su3::Triplet => 1,
        Su3::AntiTriplet => 2,
        Su3::Sextet => 3,
        Su3::AntiSextet => 4,
        Su3::Octet => 5,
    };
    let yk = (y * 3600.0).round() as i64;
    let a = (code(su3), su2, yk);
    let b = (code(su3.conjugate()), su2, -yk);
    a.min(b)
}

/// True if the active fields form supersymmetric multiplets: one gaugino per
/// gauge group, and the remaining Weyl fermions match the complex scalars
/// one-to-one up to conjugation.
pub fn is_supersymmetric(active: &[&Field]) -> bool {
    if active.iter().any(|f| f.kind == Kind::RealScalar) {
        return false;
    }
    let mut fermions: HashMap<(u8, u32, i64), i64> = HashMap::new();
    for f in active {
        for (c, w, y, n) in f.weyl_components() {
            *fermions.entry(charge_key(c, w, y)).or_default() += n as i64;
        }
    }
    let gauginos = [
        charge_key(Su3::Octet, 1, 0.0),
        charge_key(Su3::Singlet, 3, 0.0),
        charge_key(Su3::Singlet, 1, 0.0),
    ];
    for g in gauginos {
        match fermions.get_mut(&g) {
            Some(n) if *n >= 1 => *n -= 1,
            _ => return false,
        }
    }
    let mut scalars: HashMap<(u8, u32, i64), i64> = HashMap::new();
    for f in active.iter().filter(|f| f.kind == Kind::ComplexScalar) {
        *scalars.entry(charge_key(f.su3, f.su2, f.y)).or_default() += f.copies as i64;
    }
    fermions.retain(|_, n| *n != 0);
    fermions == scalars
}

/// Two-loop matrix of a supersymmetric theory, from its chiral superfields
/// (one complex scalar per chiral superfield).
fn superfield_matrix(active: &[&Field]) -> [[f64; 3]; 3] {
    let mut b = [[0.0; 3]; 3];
    for i in 0..3 {
        b[i][i] = -6.0 * GAUGE_CASIMIR[i] * GAUGE_CASIMIR[i];
    }
    for f in active.iter().filter(|f| f.kind == Kind::ComplexScalar) {
        let s2 = f.dynkin();
        let c2 = f.casimirs();
        for i in 0..3 {
            for j in 0..3 {
                let diag = if i == j { 2.0 * GAUGE_CASIMIR[i] } else { 0.0 };
                b[i][j] += (4.0 * c2[j] + diag) * s2[i];
            }
        }
    }
    b
}

/// Two-loop matrix for the active fields, choosing the supersymmetric form when it applies.
pub fn active_matrix(active: &[&Field]) -> [[f64; 3]; 3] {
    if is_supersymmetric(active) {
        superfield_matrix(active)
    } else {
        two_loop_matrix(active.iter().copied())
    }
}

#[derive(Debug, Clone)]
struct Segment {
    lo: f64,
    hi: f64,
    b: [f64; 3],
    bij: [[f64; 3]; 3],
    yukawa: Option<YukawaCoeffs>,
}

/// State along the running: ln Q, 1/alpha_i, y_t.
#[derive(Debug, Clone, Copy)]
pub struct Point {
    pub ln_q: f64,
    pub inv: [f64; 3],
    pub y_top: f64,
}

#[derive(Debug, Clone)]
pub struct TwoLoop {
    segments: Vec<Segment>,
}

impl TwoLoop {
    pub fn new(fields: &[Field], top_yukawa: TopYukawa) -> Self {
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
            .map(|w| {
                let active: Vec<&Field> =
                    fields.iter().filter(|f| f.mass <= w[0] * (1.0 + 1e-9)).collect();
                let susy = is_supersymmetric(&active);
                let yukawa = match top_yukawa {
                    TopYukawa::None => None,
                    TopYukawa::Auto if susy => Some(MSSM_YUKAWA),
                    TopYukawa::Auto => Some(SM_YUKAWA),
                };
                Segment {
                    lo: w[0],
                    hi: w[1],
                    b: beta_coefficients(active.iter().copied()),
                    bij: active_matrix(&active),
                    yukawa,
                }
            })
            .collect();
        Self { segments }
    }

    fn derivative(seg: &Segment, inv: [f64; 3], y: f64) -> ([f64; 3], f64) {
        let alpha = inv.map(|x| 1.0 / x);
        let mut d_inv = [0.0; 3];
        for i in 0..3 {
            let two_loop: f64 = (0..3).map(|j| seg.bij[i][j] * alpha[j]).sum();
            d_inv[i] = -seg.b[i] / (2.0 * PI) - two_loop / (8.0 * PI * PI);
        }
        let mut d_y = 0.0;
        if let Some(yk) = seg.yukawa {
            for i in 0..3 {
                d_inv[i] += yk.d[i] * y * y / (32.0 * PI * PI * PI);
            }
            let gauge: f64 = (0..3).map(|i| yk.c[i] * 4.0 * PI * alpha[i]).sum();
            d_y = y / (16.0 * PI * PI) * (yk.a * y * y - gauge);
        }
        (d_inv, d_y)
    }

    /// Integrate with RK4 from M_Z up to Q_MAX, stopping early if a coupling
    /// becomes non-perturbative.
    pub fn trajectory(&self, inv_mz: [f64; 3]) -> Vec<Point> {
        const STEP: f64 = 0.02;
        let mut inv = inv_mz;
        let mut y = Y_TOP_MZ;
        let mut points = vec![Point { ln_q: M_Z.ln(), inv, y_top: y }];
        for seg in &self.segments {
            let (t0, t1) = (seg.lo.ln(), seg.hi.ln());
            let n = ((t1 - t0) / STEP).ceil().max(1.0) as usize;
            let h = (t1 - t0) / n as f64;
            for step in 0..n {
                let add = |a: [f64; 3], k: [f64; 3], s: f64| [a[0] + s * k[0], a[1] + s * k[1], a[2] + s * k[2]];
                let (k1, l1) = Self::derivative(seg, inv, y);
                let (k2, l2) = Self::derivative(seg, add(inv, k1, h / 2.0), y + h / 2.0 * l1);
                let (k3, l3) = Self::derivative(seg, add(inv, k2, h / 2.0), y + h / 2.0 * l2);
                let (k4, l4) = Self::derivative(seg, add(inv, k3, h), y + h * l3);
                for i in 0..3 {
                    inv[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
                }
                y += h / 6.0 * (l1 + 2.0 * l2 + 2.0 * l3 + l4);
                if inv.iter().any(|&x| !(x > 0.0)) || !(y.abs() < 10.0) {
                    return points;
                }
                points.push(Point { ln_q: t0 + (step + 1) as f64 * h, inv, y_top: y });
            }
        }
        points
    }

    /// First point where couplings i and j meet, by linear interpolation.
    pub fn crossing(points: &[Point], i: usize, j: usize) -> Option<Point> {
        points.windows(2).find_map(|w| {
            let (d0, d1) = (w[0].inv[i] - w[0].inv[j], w[1].inv[i] - w[1].inv[j]);
            if d0 == 0.0 {
                return Some(w[0]);
            }
            if d0.signum() == d1.signum() {
                return None;
            }
            let s = d0 / (d0 - d1);
            let lerp = |a: f64, b: f64| a + s * (b - a);
            Some(Point {
                ln_q: lerp(w[0].ln_q, w[1].ln_q),
                inv: [0, 1, 2].map(|k| lerp(w[0].inv[k], w[1].inv[k])),
                y_top: lerp(w[0].y_top, w[1].y_top),
            })
        })
    }

    /// Solve for the held-out 1/alpha_k(M_Z) such that all three couplings
    /// meet where the two inputs meet. `guess` seeds the bracket (e.g. the
    /// one-loop answer). Returns (1/alpha_k(M_Z), crossing point).
    pub fn solve(
        &self,
        inv_measured: [f64; 3],
        inputs: (usize, usize),
        k: usize,
        guess: f64,
    ) -> Option<(f64, Point)> {
        let (i, j) = inputs;
        let mismatch = |x: f64| {
            let mut inv = inv_measured;
            inv[k] = x;
            let p = Self::crossing(&self.trajectory(inv), i, j)?;
            Some((p.inv[k] - p.inv[i], p))
        };

        // Expand a bracket around the guess until the mismatch changes sign.
        let mut width = 0.05 * guess.abs().max(1.0);
        let (mut lo, mut hi) = (None, None);
        for _ in 0..12 {
            let (a, b) = ((guess - width).max(0.5), guess + width);
            if let (Some((fa, _)), Some((fb, _))) = (mismatch(a), mismatch(b)) {
                if fa.signum() != fb.signum() {
                    lo = Some((a, fa));
                    hi = Some((b, fb));
                    break;
                }
            }
            width *= 2.0;
        }
        let (mut a, fa) = lo?;
        let (mut b, _) = hi?;
        for _ in 0..60 {
            let m = 0.5 * (a + b);
            let (fm, _) = mismatch(m)?;
            if fm.signum() == fa.signum() {
                a = m;
            } else {
                b = m;
            }
        }
        let x = 0.5 * (a + b);
        mismatch(x).map(|(_, p)| (x, p))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::couplings::unification::inverse_couplings_mz;
    use crate::model::field::{Kind, Su3};
    use crate::model::running::Running;
    use crate::model::sm::standard_model;

    fn assert_matrix(got: [[f64; 3]; 3], want: [[f64; 3]; 3]) {
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (got[i][j] - want[i][j]).abs() < 1e-9,
                    "b[{}][{}] = {}, want {}",
                    i + 1,
                    j + 1,
                    got[i][j],
                    want[i][j]
                );
            }
        }
    }

    #[test]
    fn test_sm_two_loop_matrix() {
        assert_matrix(
            two_loop_matrix(&standard_model()),
            [
                [199.0 / 50.0, 27.0 / 10.0, 44.0 / 5.0],
                [9.0 / 10.0, 35.0 / 6.0, 12.0],
                [11.0 / 10.0, 9.0 / 2.0, -26.0],
            ],
        );
    }

    fn mssm_fields() -> Vec<Field> {
        let mut fields = standard_model();
        fields.extend([
            Field::new("gluino", Kind::Weyl, Su3::Octet, 1, 0.0, 1),
            Field::new("wino", Kind::Weyl, Su3::Singlet, 3, 0.0, 1),
            Field::new("bino", Kind::Weyl, Su3::Singlet, 1, 0.0, 1),
            Field::new("higgsinos", Kind::Dirac, Su3::Singlet, 2, 0.5, 1),
            Field::new("H_d", Kind::ComplexScalar, Su3::Singlet, 2, -0.5, 1),
            Field::new("sQ", Kind::ComplexScalar, Su3::Triplet, 2, 1.0 / 6.0, 3),
            Field::new("su", Kind::ComplexScalar, Su3::AntiTriplet, 1, -2.0 / 3.0, 3),
            Field::new("sd", Kind::ComplexScalar, Su3::AntiTriplet, 1, 1.0 / 3.0, 3),
            Field::new("sL", Kind::ComplexScalar, Su3::Singlet, 2, -0.5, 3),
            Field::new("se", Kind::ComplexScalar, Su3::Singlet, 1, 1.0, 3),
        ]);
        fields
    }

    #[test]
    fn test_supersymmetry_detection() {
        let mssm = mssm_fields();
        assert!(is_supersymmetric(&mssm.iter().collect::<Vec<_>>()));
        let sm = standard_model();
        assert!(!is_supersymmetric(&sm.iter().collect::<Vec<_>>()));
        // Drop the selectrons: no longer supersymmetric.
        let partial: Vec<&Field> = mssm.iter().filter(|f| f.name != "se").collect();
        assert!(!is_supersymmetric(&partial));
    }

    #[test]
    fn test_mssm_two_loop_matrix() {
        let fields = mssm_fields();
        assert_matrix(
            active_matrix(&fields.iter().collect::<Vec<_>>()),
            [
                [199.0 / 25.0, 27.0 / 5.0, 88.0 / 5.0],
                [9.0 / 5.0, 25.0, 24.0],
                [11.0 / 5.0, 9.0, 14.0],
            ],
        );
    }

    #[test]
    fn test_reduces_to_one_loop() {
        // With b_ij and Yukawas switched off, RK4 must match the analytic one-loop result.
        let fields = standard_model();
        let mut two = TwoLoop::new(&fields, TopYukawa::None);
        for s in &mut two.segments {
            s.bij = [[0.0; 3]; 3];
        }
        let inv = inverse_couplings_mz();
        let last = *two.trajectory(inv).last().unwrap();
        let one = Running::new(&fields).inverse_at(inv, last.ln_q.exp());
        for i in 0..3 {
            assert!((last.inv[i] - one[i]).abs() < 1e-8, "{} vs {}", last.inv[i], one[i]);
        }
    }

    #[test]
    fn test_solution_unifies() {
        let fields = standard_model();
        let two = TwoLoop::new(&fields, TopYukawa::Auto);
        let (x, p) = two.solve(inverse_couplings_mz(), (0, 1), 2, 14.0).unwrap();
        assert!((p.inv[2] - p.inv[0]).abs() < 1e-6);
        assert!(x > 10.0 && x < 20.0, "SM two-loop 1/alpha_s = {}", x);
    }
}
