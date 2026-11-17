//! Brute-force search for non-supersymmetric amendments that unify.
//!
//! Every combination of up to `max_fields` multiplets from a fixed menu, each at
//! a mass from a grid, is screened with the analytic one-loop prediction. The
//! survivors are re-evaluated at two loops and kept if they predict alpha_s,
//! survive the Super-K proton limit, unify below 1e18 GeV and stay perturbative.
//! Menu entries are vectorlike fermions, Weyl fermions in real representations,
//! or scalars, so every combination is anomaly-free by construction.

use std::fmt::Write as _;

use std::sync::atomic::{AtomicUsize, Ordering};

use rayon::prelude::*;

use super::field::{Field, Kind, Su3};
use super::predict::{predict, Prediction};
use super::proton::{self, ProtonDecay};
use super::running::Running;
use super::spec::{Amendment, Base, TopYukawa};
use super::two_loop::TwoLoop;
use crate::couplings::unification::inverse_couplings_mz;

/// A candidate multiplet.
#[derive(Debug, Clone, Copy)]
pub struct MenuItem {
    pub label: &'static str,
    pub kind: Kind,
    pub su3: Su3,
    pub su2: u32,
    pub y: f64,
    /// Where this multiplet appears in the literature.
    pub note: &'static str,
}

impl MenuItem {
    pub fn field(&self, copies: u32, mass: f64) -> Field {
        Field::new(self.label, self.kind, self.su3, self.su2, self.y, copies).with_mass(mass)
    }

    /// Number of new particle states (complex components count once, Dirac twice).
    pub fn states(&self) -> u32 {
        let d = self.su3.dim() as u32 * self.su2;
        match self.kind {
            Kind::Dirac => 2 * d,
            _ => d,
        }
    }
}

pub fn menu() -> Vec<MenuItem> {
    use Kind::*;
    use Su3::*;
    let m = |label, kind, su3, su2, y, note| MenuItem { label, kind, su3, su2, y, note };
    vec![
        m("Weyl (1,3,0)", Weyl, Singlet, 3, 0.0, "wino-like triplet"),
        m("Weyl (8,1,0)", Weyl, Octet, 1, 0.0, "gluino-like octet"),
        m("Weyl (1,5,0)", Weyl, Singlet, 5, 0.0, "minimal dark matter quintuplet"),
        m("Dirac (1,2,1/2)", Dirac, Singlet, 2, 0.5, "higgsino-like doublet"),
        m("Dirac (3,2,1/6)", Dirac, Triplet, 2, 1.0 / 6.0, "vectorlike quark doublet"),
        m("Dirac (3,1,2/3)", Dirac, Triplet, 1, 2.0 / 3.0, "vectorlike top partner"),
        m("Dirac (3,1,-1/3)", Dirac, Triplet, 1, -1.0 / 3.0, "vectorlike bottom partner"),
        m("Dirac (1,1,1)", Dirac, Singlet, 1, 1.0, "vectorlike charged lepton"),
        m("Dirac (1,3,1)", Dirac, Singlet, 3, 1.0, "vectorlike triplet lepton"),
        m("Dirac (3,3,-1/3)", Dirac, Triplet, 3, -1.0 / 3.0, "vectorlike quark triplet"),
        m("Scalar (1,2,1/2)", ComplexScalar, Singlet, 2, 0.5, "second Higgs doublet"),
        m("Scalar (3,1,-1/3)", ComplexScalar, Triplet, 1, -1.0 / 3.0, "leptoquark S1"),
        m("Scalar (3,2,1/6)", ComplexScalar, Triplet, 2, 1.0 / 6.0, "leptoquark R2-tilde"),
        m("Scalar (3,2,7/6)", ComplexScalar, Triplet, 2, 7.0 / 6.0, "leptoquark R2"),
        m("Scalar (3,3,-1/3)", ComplexScalar, Triplet, 3, -1.0 / 3.0, "leptoquark S3"),
        m("Scalar (8,2,1/2)", ComplexScalar, Octet, 2, 0.5, "Manohar-Wise colour octet"),
        m("Scalar (1,3,1)", ComplexScalar, Singlet, 3, 1.0, "type-II seesaw triplet"),
        m("Scalar (6,1,1/3)", ComplexScalar, Sextet, 1, 1.0 / 3.0, "diquark"),
        m("Real scalar (1,3,0)", RealScalar, Singlet, 3, 0.0, "real triplet"),
        m("Real scalar (8,1,0)", RealScalar, Octet, 1, 0.0, "real colour octet"),
    ]
}

/// Default mass grid (GeV): the lightest is about where the LHC stops.
pub const DEFAULT_MASSES: [f64; 7] = [1e3, 1e4, 1e5, 1e6, 1e8, 1e10, 1e12];

/// Unification above this scale would be in the quantum-gravity regime.
pub const MAX_M_GUT: f64 = 1e18;

#[derive(Debug, Clone)]
pub struct SearchConfig {
    pub max_fields: usize,
    pub masses: Vec<f64>,
    pub loops: u8,
    /// One-loop screen: keep candidates whose alpha_s is within this relative distance.
    pub screen: f64,
    /// Accept candidates with |pull| at or below this.
    pub accept_sigma: f64,
}

impl Default for SearchConfig {
    fn default() -> Self {
        Self {
            max_fields: 3,
            masses: DEFAULT_MASSES.to_vec(),
            loops: 2,
            screen: 0.15,
            accept_sigma: 1.0,
        }
    }
}

/// An accepted candidate.
#[derive(Debug, Clone)]
pub struct Hit {
    /// (menu index, mass) for each added multiplet, sorted.
    pub picks: Vec<(usize, f64)>,
    pub amendment: Amendment,
    pub prediction: Prediction,
    pub proton: ProtonDecay,
    pub new_states: u32,
    pub mass_scales: usize,
    pub lightest_coloured: Option<f64>,
}

impl Hit {
    pub fn fields(&self) -> usize {
        self.picks.len()
    }
}

#[derive(Debug, Clone)]
pub struct SearchResult {
    pub evaluated: usize,
    pub screened: usize,
    pub hits: Vec<Hit>,
    /// Two-loop prediction for the SM alone, for reference.
    pub baseline: Option<Prediction>,
}

/// Build an amendment from picks, merging identical (item, mass) picks into copies.
pub fn amendment_for(menu: &[MenuItem], picks: &[(usize, f64)]) -> Amendment {
    let mut added: Vec<Field> = Vec::new();
    let mut names = Vec::new();
    let mut k = 0;
    while k < picks.len() {
        let mut n = 1;
        while k + n < picks.len() && picks[k + n] == picks[k] {
            n += 1;
        }
        let (item, mass) = picks[k];
        added.push(menu[item].field(n as u32, mass));
        names.push(format!(
            "{}{} @ {:.0e} GeV",
            if n > 1 { format!("{} x ", n) } else { String::new() },
            menu[item].label,
            mass
        ));
        k += n;
    }
    Amendment {
        name: format!("SM + {}", names.join(" + ")),
        description: Some("Found by `reality search`.".to_string()),
        base: Base::Sm,
        inputs: (0, 1),
        predicted: 2,
        top_yukawa: TopYukawa::Auto,
        tan_beta: None,
        added,
        heavy: Vec::new(),
    }
}

/// Every multiset of `size` items from 0..n whose first (smallest) element is
/// `first`, in non-decreasing order. Splitting on `first` lets the search
/// spread the enumeration over threads without storing every candidate.
fn multisets_from(first: usize, n: usize, size: usize, mut visit: impl FnMut(&[usize])) {
    if size == 0 || first >= n {
        return;
    }
    let mut idx = vec![first; size];
    if size == 1 {
        visit(&idx);
        return;
    }
    multisets(n - first, size - 1, |rest| {
        for (k, &r) in rest.iter().enumerate() {
            idx[k + 1] = first + r;
        }
        visit(&idx);
    });
}

/// Every multiset of `size` items from 0..n, in non-decreasing order.
fn multisets(n: usize, size: usize, mut visit: impl FnMut(&[usize])) {
    let mut idx = vec![0usize; size];
    if size == 0 || n == 0 {
        return;
    }
    loop {
        visit(&idx);
        let mut pos = size;
        while pos > 0 && idx[pos - 1] == n - 1 {
            pos -= 1;
        }
        if pos == 0 {
            return;
        }
        idx[pos - 1] += 1;
        for p in pos..size {
            idx[p] = idx[pos - 1];
        }
    }
}

fn full_prediction(amendment: &Amendment, loops: u8) -> Option<Prediction> {
    let fields = amendment.all_fields();
    let running = Running::new(&fields);
    let inv = inverse_couplings_mz();
    let one = predict(&running, amendment.inputs, amendment.predicted, [0.0; 3])?;
    let p = if loops >= 2 {
        let two = TwoLoop::new(&fields, amendment.top_yukawa);
        let (x, pt) = two.solve(inv, amendment.inputs, amendment.predicted, one.inverse_predicted, [0.0; 3])?;
        Prediction::new(2, amendment.inputs, amendment.predicted, pt.ln_q.exp(), pt.inv[0], x)
    } else {
        one
    };
    match running.landau_pole(inv) {
        Some(pole) if pole <= p.m_gut => None,
        _ => Some(p),
    }
}

pub fn search(cfg: &SearchConfig) -> SearchResult {
    let menu = menu();
    let items: Vec<(usize, f64)> = (0..menu.len())
        .flat_map(|i| cfg.masses.iter().map(move |&m| (i, m)))
        .collect();
    let inv = inverse_couplings_mz();
    let alpha_measured = 1.0 / inv[2];

    // Stage 1: analytic one-loop screen. Candidates are generated inside each
    // thread (split by their first item), so they are never all held in memory.
    let evaluated = AtomicUsize::new(0);
    let item_fields: Vec<Field> = items.iter().map(|&(i, m)| menu[i].field(1, m)).collect();
    let sm = super::sm::standard_model();
    let mut survivors: Vec<Vec<(usize, f64)>> = Vec::new();
    for size in 1..=cfg.max_fields {
        let batch: Vec<Vec<(usize, f64)>> = (0..items.len())
            .into_par_iter()
            .flat_map_iter(|first| {
                let mut kept = Vec::new();
                let mut count = 0;
                let mut fields = sm.clone();
                multisets_from(first, items.len(), size, |idx| {
                    count += 1;
                    fields.truncate(sm.len());
                    fields.extend(idx.iter().map(|&k| item_fields[k].clone()));
                    let running = Running::new(&fields);
                    let pass = predict(&running, (0, 1), 2, [0.0; 3]).map_or(false, |p| {
                        (p.alpha_predicted() / alpha_measured - 1.0).abs() <= cfg.screen
                            && p.m_gut <= MAX_M_GUT
                    });
                    if pass {
                        kept.push(idx.iter().map(|&k| items[k]).collect());
                    }
                });
                evaluated.fetch_add(count, Ordering::Relaxed);
                kept
            })
            .collect();
        survivors.extend(batch);
    }
    let evaluated = evaluated.into_inner();
    let screened = survivors.len();

    // Stage 2: full prediction, proton decay and acceptance, in parallel.
    // Collecting an indexed parallel iterator keeps the input order, so the
    // result is identical to a sequential run.
    let mut hits: Vec<Hit> = survivors
        .into_par_iter()
        .filter_map(|picks| {
            let amendment = amendment_for(&menu, &picks);
            let prediction = full_prediction(&amendment, cfg.loops)?;
            if prediction.pull.abs() > cfg.accept_sigma || prediction.m_gut > MAX_M_GUT {
                return None;
            }
            let proton = ProtonDecay::estimate(prediction.m_gut, prediction.inverse_alpha_gut);
            if proton.lifetime_years < proton::SUPER_K_LIMIT_YEARS {
                return None;
            }
            let mut scales: Vec<f64> = picks.iter().map(|p| p.1).collect();
            scales.dedup();
            let lightest_coloured = picks
                .iter()
                .filter(|&&(i, _)| menu[i].su3 != Su3::Singlet)
                .map(|p| p.1)
                .reduce(f64::min);
            Some(Hit {
                new_states: picks.iter().map(|&(i, _)| menu[i].states()).sum(),
                mass_scales: scales.len(),
                lightest_coloured,
                picks,
                amendment,
                prediction,
                proton,
            })
        })
        .collect();

    hits.sort_by(|a, b| {
        (a.fields(), a.mass_scales, a.new_states)
            .cmp(&(b.fields(), b.mass_scales, b.new_states))
            .then(a.prediction.pull.abs().partial_cmp(&b.prediction.pull.abs()).unwrap())
    });

    let sm = Amendment {
        added: Vec::new(),
        ..amendment_for(&menu, &[])
    };
    SearchResult {
        evaluated,
        screened,
        hits,
        baseline: full_prediction(&sm, cfg.loops),
    }
}

/// Write a rational hypercharge such as -1/3 or 7/6.
fn rational(y: f64) -> String {
    for d in 1..=36 {
        let n = (y * d as f64).round();
        if (y * d as f64 - n).abs() < 1e-9 {
            return if d == 1 { format!("{}", n as i64) } else { format!("{}/{}", n as i64, d) };
        }
    }
    format!("{}", y)
}

fn kind_name(kind: Kind) -> &'static str {
    match kind {
        Kind::Weyl => "weyl",
        Kind::Dirac => "dirac",
        Kind::ComplexScalar => "complex_scalar",
        Kind::RealScalar => "real_scalar",
        Kind::Chiral => "chiral",
        Kind::Vector => "vector",
        Kind::VectorSuperfield => "vector_superfield",
    }
}

/// The amendment as a model file that `reality amend` can read.
pub fn to_toml(hit: &Hit) -> String {
    let a = &hit.amendment;
    let mut s = String::new();
    let _ = writeln!(s, "# Found by `reality search`: alpha_s(M_Z) = {:.4} at two loops, pull {:+.2},", hit.prediction.alpha_predicted(), hit.prediction.pull);
    let _ = writeln!(s, "# M_GUT = {:.2e} GeV, proton lifetime ~{:.1e} years.", hit.prediction.m_gut, hit.proton.lifetime_years);
    let _ = writeln!(s, "name = {:?}", a.name);
    let _ = writeln!(s, "inputs = [\"alpha_1\", \"alpha_2\"]");
    for f in &a.added {
        let _ = writeln!(s, "\n[[field]]");
        let _ = writeln!(s, "name = {:?}", f.name);
        let _ = writeln!(s, "kind = \"{}\"", kind_name(f.kind));
        let _ = writeln!(s, "su3 = \"{}\"", f.su3.label());
        let _ = writeln!(s, "su2 = {}", f.su2);
        let _ = writeln!(s, "y = \"{}\"", rational(f.y));
        if f.copies > 1 {
            let _ = writeln!(s, "copies = {}", f.copies);
        }
        let _ = writeln!(s, "mass = {:e}", f.mass);
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::gates::{check, Verdict};
    use crate::model::sm::standard_model;

    #[test]
    fn test_menu_items_are_consistent() {
        // Each entry must pass every gate on its own; gates are additive,
        // so then every combination passes too.
        for item in menu() {
            let mut fields = standard_model();
            fields.push(item.field(1, 1000.0));
            for g in check(&fields) {
                assert_eq!(g.verdict, Verdict::Pass, "{}: {} {}", item.label, g.name, g.detail);
            }
        }
    }

    #[test]
    fn test_multisets_count() {
        let mut n = 0;
        multisets(5, 3, |_| n += 1);
        assert_eq!(n, 35); // C(5 + 3 - 1, 3)
    }

    #[test]
    fn test_multisets_from_partitions_multisets() {
        // Splitting by first element visits exactly the same multisets, in the same order.
        let mut all = Vec::new();
        multisets(6, 3, |m| all.push(m.to_vec()));
        let mut split = Vec::new();
        for first in 0..6 {
            multisets_from(first, 6, 3, |m| split.push(m.to_vec()));
        }
        assert_eq!(all, split);
    }

    #[test]
    fn test_toml_round_trip() {
        let menu = menu();
        let picks = vec![(0, 1e3), (0, 1e3), (13, 1e6)];
        let a = amendment_for(&menu, &picks);
        assert_eq!(a.added.len(), 2);
        assert_eq!(a.added[0].copies, 2);
        let hit = Hit {
            prediction: full_prediction(&a, 1).unwrap(),
            proton: ProtonDecay::estimate(1e16, 25.0),
            new_states: 0,
            mass_scales: 2,
            lightest_coloured: None,
            picks,
            amendment: a.clone(),
        };
        let back = Amendment::from_toml(&to_toml(&hit)).unwrap();
        assert_eq!(back.added.len(), 2);
        for (x, y) in back.added.iter().zip(&a.added) {
            assert_eq!((x.su3, x.su2, x.copies), (y.su3, y.su2, y.copies));
            assert!((x.y - y.y).abs() < 1e-12 && (x.mass / y.mass - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn test_small_search_accepts_only_good_candidates() {
        // Everything at 1 TeV: the SM baseline fails, and every accepted
        // candidate meets the alpha_s and proton criteria.
        let cfg = SearchConfig {
            max_fields: 3,
            masses: vec![1e3],
            ..SearchConfig::default()
        };
        let result = search(&cfg);
        assert!(result.evaluated > 0);
        let baseline = result.baseline.unwrap();
        assert!(baseline.pull.abs() > 5.0, "SM should fail: pull {}", baseline.pull);
        assert!(result.hits.iter().all(|h| h.prediction.pull.abs() <= 1.0));
        assert!(result
            .hits
            .iter()
            .all(|h| h.proton.lifetime_years >= proton::SUPER_K_LIMIT_YEARS));
    }
}
