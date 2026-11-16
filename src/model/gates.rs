use super::field::{su2_index, Field, Su3};

const TOL: f64 = 1e-9;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verdict {
    Pass,
    /// Allowed in principle, but strongly constrained by experiment.
    Warn,
    Fail,
}

/// A consistency check that needs no experimental data.
#[derive(Debug, Clone)]
pub struct Gate {
    pub name: &'static str,
    pub verdict: Verdict,
    pub detail: String,
}

fn sum_gate(name: &'static str, value: f64) -> Gate {
    let ok = value.abs() < TOL;
    Gate {
        name,
        verdict: if ok { Verdict::Pass } else { Verdict::Fail },
        detail: format!("sum = {}", if ok { "0".to_string() } else { format!("{:.6}", value) }),
    }
}

/// Gauge anomaly cancellation, the Witten SU(2) anomaly, and integer charges
/// for colour-singlet fermions.
pub fn check(fields: &[Field]) -> Vec<Gate> {
    let weyl: Vec<(Su3, u32, f64, f64)> = fields
        .iter()
        .flat_map(|f| f.weyl_components())
        .map(|(c, w, y, n)| (c, w, y, n as f64))
        .collect();

    let total = |g: &dyn Fn(Su3, u32, f64) -> f64| -> f64 {
        weyl.iter().map(|&(c, w, y, n)| n * g(c, w, y)).sum()
    };

    let mut gates = vec![
        sum_gate("U(1)_Y^3", total(&|c, w, y| c.dim() * w as f64 * y.powi(3))),
        sum_gate("gravity^2 x U(1)_Y", total(&|c, w, y| c.dim() * w as f64 * y)),
        sum_gate("SU(3)^2 x U(1)_Y", total(&|c, w, y| c.index() * w as f64 * y)),
        sum_gate("SU(2)^2 x U(1)_Y", total(&|c, w, y| su2_index(w) * c.dim() * y)),
        sum_gate("SU(3)^3", total(&|c, w, _| c.cubic_anomaly() * w as f64)),
    ];

    // Witten anomaly: the number of SU(2) multiplets with odd 2T(R),
    // counted with colour multiplicity, must be even.
    let witten: f64 = total(&|c, w, _| c.dim() * (2.0 * su2_index(w)));
    let even = (witten.round() as i64) % 2 == 0;
    gates.push(Gate {
        name: "Witten SU(2)",
        verdict: if even { Verdict::Pass } else { Verdict::Fail },
        detail: format!("{} doublet-equivalents ({})", witten.round(), if even { "even" } else { "odd" }),
    });

    // Colour-singlet fermions with fractional electric charge would be stable
    // relics; searches constrain them heavily.
    let mut fractional = Vec::new();
    for f in fields {
        if f.su3 != Su3::Singlet || f.weyl_components().is_empty() {
            continue;
        }
        let j2 = f.su2 as i32 - 1; // 2j
        for m2 in (-j2..=j2).step_by(2) {
            let q = 0.5 * m2 as f64 + f.y;
            if (q - q.round()).abs() > TOL {
                fractional.push(format!("{} (Q = {:.3})", f.name, q));
                break;
            }
        }
    }
    gates.push(Gate {
        name: "integer charges",
        verdict: if fractional.is_empty() { Verdict::Pass } else { Verdict::Warn },
        detail: if fractional.is_empty() {
            "all colour-singlet fermions".to_string()
        } else {
            format!("fractional: {}", fractional.join(", "))
        },
    });

    gates
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::field::Kind;
    use crate::model::sm::standard_model;

    #[test]
    fn test_sm_passes_all() {
        for g in check(&standard_model()) {
            assert_eq!(g.verdict, Verdict::Pass, "{}: {}", g.name, g.detail);
        }
    }

    #[test]
    fn test_single_higgsino_fails() {
        let mut fields = standard_model();
        fields.push(Field::new("higgsino", Kind::Weyl, Su3::Singlet, 2, 0.5, 1));
        let gates = check(&fields);
        let failed: Vec<_> = gates.iter().filter(|g| g.verdict == Verdict::Fail).map(|g| g.name).collect();
        assert!(failed.contains(&"U(1)_Y^3"));
        assert!(failed.contains(&"Witten SU(2)"));
    }

    #[test]
    fn test_vectorlike_is_anomaly_free() {
        let mut fields = standard_model();
        fields.push(Field::new("vlq", Kind::Dirac, Su3::Triplet, 2, 1.0 / 6.0, 1));
        assert!(check(&fields).iter().all(|g| g.verdict == Verdict::Pass));
    }

    #[test]
    fn test_missing_quark_breaks_su3_cubed() {
        let mut fields = standard_model();
        fields.retain(|f| f.name != "d^c");
        let gates = check(&fields);
        assert!(gates.iter().any(|g| g.name == "SU(3)^3" && g.verdict == Verdict::Fail));
    }

    #[test]
    fn test_fractional_charge_warns() {
        let mut fields = standard_model();
        fields.push(Field::new("odd", Kind::Dirac, Su3::Singlet, 1, 0.5, 1));
        let gates = check(&fields);
        assert!(gates.iter().any(|g| g.name == "integer charges" && g.verdict == Verdict::Warn));
    }
}
