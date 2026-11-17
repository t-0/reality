//! Robustness: how much of the nearby model space still works.
//!
//! A single best-fit point says little about how tuned a model is. Here each
//! distinct mass scale of the added fields (fields sharing a mass move together)
//! is drawn log-uniformly within `spread_decades` of its value, and
//! the GUT-scale fields are drawn log-uniformly within their ranges. The
//! robustness is the fraction of draws that still pass: no new particle below
//! `min_mass` (roughly where colliders stop), |pull| <= accept_sigma, the proton
//! outlives the Super-K limit, M_GUT <= 1e18 GeV, and all couplings stay
//! perturbative. Draws below `min_mass` count as failures rather than being
//! clamped, so a model at the edge of the excluded region is not favoured.

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;

use super::field::{Field, Kind, Su3};
use super::predict::evaluate;
use super::proton;
use super::search::MAX_M_GUT;
use super::spec::{Amendment, HeavyField};

#[derive(Debug, Clone)]
pub struct RobustConfig {
    pub samples: usize,
    /// Each added mass moves log-uniformly within this many decades either way.
    pub spread_decades: f64,
    /// Draws with any added mass below this (GeV) fail as excluded by colliders.
    pub min_mass: f64,
    /// Generic GUT-field ratios are drawn from [1/threshold_range, threshold_range]
    /// when the model declares no [[heavy]] fields. 1.0 switches them off.
    pub threshold_range: f64,
    pub accept_sigma: f64,
    pub loops: u8,
    pub seed: u64,
}

impl Default for RobustConfig {
    fn default() -> Self {
        Self {
            samples: 200,
            spread_decades: 1.0,
            min_mass: 1000.0,
            threshold_range: 10.0,
            accept_sigma: 1.0,
            loops: 2,
            seed: 1,
        }
    }
}

#[derive(Debug, Clone)]
pub struct RobustReport {
    pub samples: usize,
    /// Fraction of draws that pass everything.
    pub pass: f64,
    /// Statistical (binomial) 1-sigma error on `pass`.
    pub pass_error: f64,
    /// Fraction of draws with a new particle below `min_mass`.
    pub fail_collider: f64,
    /// Fractions failing each criterion among draws above `min_mass`
    /// (a draw can fail several).
    pub fail_alpha: f64,
    pub fail_proton: f64,
    /// Draws with no prediction at all (no unification, Landau pole, M_GUT too high).
    pub fail_unify: f64,
    /// 5th, 50th and 95th percentiles of the predicted alpha_s over draws.
    pub alpha_quantiles: Option<[f64; 3]>,
    /// True if generic SU(5) GUT fields were added.
    pub generic_thresholds: bool,
    /// True if the model is supersymmetric at high energy (SUSY SU(5) fields used).
    pub supersymmetric: bool,
}

/// GUT-scale fields of minimal supersymmetric SU(5): the X, Y vector superfield,
/// the colour-triplet Higgs pair, and the adjoint octet and triplet, which share
/// one mass in the minimal model.
pub fn generic_susy_su5_heavy(range: f64) -> Vec<HeavyField> {
    let h = |name: &str, kind, su3, su2, y: f64, copies, ratio: (f64, f64), group: Option<&str>| HeavyField {
        field: Field::new(name, kind, su3, su2, y, copies).with_mass(0.0),
        ratio,
        group: group.map(str::to_string),
    };
    let r = (1.0 / range, range);
    vec![
        h("X, Y vector superfield", Kind::VectorSuperfield, Su3::Triplet, 2, -5.0 / 6.0, 1, (1.0, 1.0), None),
        h("colour-triplet Higgs pair", Kind::Chiral, Su3::Triplet, 1, -1.0 / 3.0, 2, r, None),
        h("adjoint octet", Kind::Chiral, Su3::Octet, 1, 0.0, 1, r, Some("sigma")),
        h("adjoint triplet", Kind::Chiral, Su3::Singlet, 3, 0.0, 1, r, Some("sigma")),
    ]
}

/// GUT-scale fields of minimal non-supersymmetric SU(5): X, Y gauge bosons
/// (defining M_X), the colour-triplet partner of the Higgs, and the colour-octet
/// and weak-triplet remnants of the adjoint Higgs.
pub fn generic_su5_heavy(range: f64) -> Vec<HeavyField> {
    let h = |name: &str, kind, su3, su2, y: f64, ratio: (f64, f64)| HeavyField {
        field: Field::new(name, kind, su3, su2, y, 1).with_mass(0.0),
        ratio,
        group: None,
    };
    let r = (1.0 / range, range);
    vec![
        h("X, Y gauge bosons", Kind::Vector, Su3::Triplet, 2, -5.0 / 6.0, (1.0, 1.0)),
        h("colour-triplet Higgs", Kind::ComplexScalar, Su3::Triplet, 1, -1.0 / 3.0, r),
        h("adjoint octet", Kind::RealScalar, Su3::Octet, 1, 0.0, r),
        h("adjoint triplet", Kind::RealScalar, Su3::Singlet, 3, 0.0, r),
    ]
}

fn log_uniform(rng: &mut StdRng, lo: f64, hi: f64) -> f64 {
    if hi <= lo {
        return lo;
    }
    (lo.ln() + rng.random::<f64>() * (hi.ln() - lo.ln())).exp()
}

pub fn robustness(amendment: &Amendment, cfg: &RobustConfig) -> RobustReport {
    let generic = amendment.heavy.is_empty() && cfg.threshold_range > 1.0;
    let all = amendment.all_fields();
    let supersymmetric = super::two_loop::is_supersymmetric(&all.iter().collect::<Vec<_>>());
    let heavy = match (generic, supersymmetric) {
        (true, true) => generic_susy_su5_heavy(cfg.threshold_range),
        (true, false) => generic_su5_heavy(cfg.threshold_range),
        (false, _) => amendment.heavy.clone(),
    };
    let spread = 10f64.powf(cfg.spread_decades);

    // Each draw: (above collider floor, passes alpha, passes proton, unified, alpha_s).
    let draws: Vec<(bool, bool, bool, bool, Option<f64>)> = (0..cfg.samples)
        .into_par_iter()
        .map(|k| {
            let mut rng = StdRng::seed_from_u64(cfg.seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ k as u64);
            // One random factor per distinct mass scale: fields that share a
            // mass are one parameter and move together.
            let mut a = amendment.clone();
            let mut scales: Vec<(f64, f64)> = Vec::new();
            for f in &mut a.added {
                let factor = match scales.iter().find(|(m, _)| (*m / f.mass - 1.0).abs() < 1e-9) {
                    Some(&(_, factor)) => factor,
                    None => {
                        let factor = log_uniform(&mut rng, 1.0 / spread, spread);
                        scales.push((f.mass, factor));
                        factor
                    }
                };
                f.mass *= factor;
            }
            let above_floor = a.added.iter().all(|f| f.mass >= cfg.min_mass * (1.0 - 1e-9));
            // One draw per heavy field, or per group for grouped fields.
            let mut groups: Vec<(String, f64)> = Vec::new();
            a.heavy = heavy
                .iter()
                .map(|h| {
                    let shared = h.group.as_ref().and_then(|g| groups.iter().find(|(n, _)| n == g).map(|x| x.1));
                    let r = shared.unwrap_or_else(|| {
                        let r = log_uniform(&mut rng, h.ratio.0, h.ratio.1);
                        if let Some(g) = &h.group {
                            groups.push((g.clone(), r));
                        }
                        r
                    });
                    HeavyField { ratio: (r, r), ..h.clone() }
                })
                .collect();
            let report = evaluate(&a, cfg.loops);
            match (report.prediction, report.proton) {
                (Some(p), Some(pd)) if p.m_gut <= MAX_M_GUT => (
                    above_floor,
                    p.pull.abs() <= cfg.accept_sigma,
                    pd.lifetime_years >= proton::SUPER_K_LIMIT_YEARS,
                    true,
                    Some(p.alpha_predicted()),
                ),
                _ => (above_floor, false, false, false, None),
            }
        })
        .collect();

    type Draw = (bool, bool, bool, bool, Option<f64>);
    let n = draws.len().max(1) as f64;
    let frac = |f: &dyn Fn(&Draw) -> bool| draws.iter().filter(|d| f(d)).count() as f64 / n;
    let mut alphas: Vec<f64> = draws.iter().filter_map(|d| d.4).collect();
    alphas.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let q = |p: f64| alphas[((alphas.len() - 1) as f64 * p).round() as usize];
    let pass = frac(&|d| d.0 && d.1 && d.2 && d.3);
    RobustReport {
        samples: draws.len(),
        pass,
        pass_error: (pass * (1.0 - pass) / n).sqrt(),
        fail_collider: frac(&|d| !d.0),
        fail_alpha: frac(&|d| d.0 && d.3 && !d.1),
        fail_proton: frac(&|d| d.0 && d.3 && !d.2),
        fail_unify: frac(&|d| d.0 && !d.3),
        alpha_quantiles: (!alphas.is_empty()).then(|| [q(0.05), q(0.5), q(0.95)]),
        generic_thresholds: generic,
        supersymmetric,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quick() -> RobustConfig {
        RobustConfig { samples: 24, ..RobustConfig::default() }
    }

    #[test]
    fn test_generic_su5_xy_matches_nonsusy_values() {
        // Non-SUSY X, Y bosons: b = -7 T(R) per (3,2,-5/6) + conjugate = (-35/2, -21/2, -7).
        let b = generic_su5_heavy(10.0)[0].field.beta_contribution();
        assert!((b[0] + 17.5).abs() < 1e-12 && (b[1] + 10.5).abs() < 1e-12 && (b[2] + 7.0).abs() < 1e-12);
    }

    #[test]
    fn test_generic_susy_su5_matches_explicit_model() {
        // The generic SUSY set reproduces models/mssm_su5.toml's heavy content.
        let explicit = Amendment::from_toml(include_str!("../../models/mssm_su5.toml")).unwrap();
        let generic = generic_susy_su5_heavy(10.0);
        assert_eq!(explicit.heavy.len(), generic.len());
        for (e, g) in explicit.heavy.iter().zip(&generic) {
            assert_eq!(e.field.beta_contribution(), g.field.beta_contribution());
            assert_eq!(e.ratio, g.ratio);
        }
    }

    #[test]
    fn test_sm_is_never_robust() {
        let sm = Amendment::from_toml("name = \"SM\"").unwrap();
        let r = robustness(&sm, &quick());
        assert_eq!(r.pass, 0.0);
    }

    #[test]
    fn test_deterministic_and_bounded() {
        let a = Amendment::from_toml(include_str!("../../models/found/quark_doublet_octet_30tev.toml")).unwrap();
        let r1 = robustness(&a, &quick());
        let r2 = robustness(&a, &quick());
        assert_eq!(r1.pass, r2.pass);
        assert!(r1.pass > 0.0 && r1.pass <= 1.0);
        let [lo, mid, hi] = r1.alpha_quantiles.unwrap();
        assert!(lo <= mid && mid <= hi);
    }

    #[test]
    fn test_floor_counts_as_failure() {
        // Split SUSY has two scales; the 1 TeV one falls below the floor about
        // half the time. All ten superpartners at that scale move together.
        let a = Amendment::from_toml(include_str!("../../models/split_susy.toml")).unwrap();
        let cfg = RobustConfig { samples: 200, threshold_range: 1.0, ..RobustConfig::default() };
        let r = robustness(&a, &cfg);
        assert!(r.fail_collider > 0.35 && r.fail_collider < 0.65, "{}", r.fail_collider);
    }

    #[test]
    fn test_no_spread_no_thresholds_reproduces_point() {
        // With nothing varied, every draw is the original model, which passes.
        let a = Amendment::from_toml(include_str!("../../models/found/quark_doublet_octet_30tev.toml")).unwrap();
        let cfg = RobustConfig { samples: 4, spread_decades: 0.0, threshold_range: 1.0, ..RobustConfig::default() };
        assert_eq!(robustness(&a, &cfg).pass, 1.0);
    }
}
