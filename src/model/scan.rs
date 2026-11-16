//! Scan an amendment over the mass of its new fields.
//!
//! All added masses are scaled together so the lightest equals the scan value,
//! which keeps any hierarchy in the model (e.g. split SUSY) intact.

use super::predict::{evaluate, Prediction};
use super::proton::ProtonDecay;
use super::spec::Amendment;

#[derive(Debug, Clone)]
pub struct ScanPoint {
    /// Mass of the lightest added field (GeV).
    pub mass: f64,
    pub prediction: Option<Prediction>,
    pub proton: Option<ProtonDecay>,
}

/// `n` masses spaced evenly in log between `from` and `to` (inclusive).
pub fn log_grid(from: f64, to: f64, n: usize) -> Vec<f64> {
    if n < 2 {
        return vec![from];
    }
    let (a, b) = (from.ln(), to.ln());
    (0..n)
        .map(|k| (a + (b - a) * k as f64 / (n - 1) as f64).exp())
        .collect()
}

impl Amendment {
    /// Copy with every added mass multiplied so the lightest becomes `lightest`.
    pub fn with_lightest_mass(&self, lightest: f64) -> Self {
        let mut out = self.clone();
        if let Some(min) = self.added.iter().map(|f| f.mass).reduce(f64::min) {
            let factor = lightest / min;
            for f in &mut out.added {
                f.mass *= factor;
            }
        }
        out
    }
}

pub fn scan(amendment: &Amendment, loops: u8, masses: &[f64]) -> Vec<ScanPoint> {
    masses
        .iter()
        .map(|&mass| {
            let report = evaluate(&amendment.with_lightest_mass(mass), loops);
            ScanPoint {
                mass,
                prediction: report.prediction,
                proton: report.proton,
            }
        })
        .collect()
}

/// Interpolate in log(mass) between two points for where `f` reaches `target`.
fn crossing_mass(m0: f64, f0: f64, m1: f64, f1: f64, target: f64) -> f64 {
    let s = (target - f0) / (f1 - f0);
    (m0.ln() + s * (m1.ln() - m0.ln())).exp()
}

/// Mass where the pull crosses zero (the best fit). If it never crosses,
/// the grid point with the smallest |pull|.
pub fn best_fit(points: &[ScanPoint]) -> Option<(f64, f64)> {
    let pulls: Vec<(f64, f64)> = points
        .iter()
        .filter_map(|p| p.prediction.as_ref().map(|pr| (p.mass, pr.pull)))
        .collect();
    for w in pulls.windows(2) {
        let ((m0, p0), (m1, p1)) = (w[0], w[1]);
        if p0 == 0.0 {
            return Some((m0, 0.0));
        }
        if p0.signum() != p1.signum() {
            return Some((crossing_mass(m0, p0, m1, p1, 0.0), 0.0));
        }
    }
    pulls
        .into_iter()
        .min_by(|a, b| a.1.abs().partial_cmp(&b.1.abs()).unwrap())
}

/// Mass ranges where |pull| <= n_sigma, with interpolated edges.
pub fn allowed_ranges(points: &[ScanPoint], n_sigma: f64) -> Vec<(f64, f64)> {
    let pulls: Vec<(f64, f64)> = points
        .iter()
        .filter_map(|p| p.prediction.as_ref().map(|pr| (p.mass, pr.pull.abs())))
        .collect();
    let mut ranges = Vec::new();
    let mut start: Option<f64> = None;
    for (k, &(m, p)) in pulls.iter().enumerate() {
        let inside = p <= n_sigma;
        match (inside, start) {
            (true, None) => {
                start = Some(if k == 0 {
                    m
                } else {
                    let (m0, p0) = pulls[k - 1];
                    crossing_mass(m0, p0, m, p, n_sigma)
                });
            }
            (false, Some(s)) => {
                let (m0, p0) = pulls[k - 1];
                ranges.push((s, crossing_mass(m0, p0, m, p, n_sigma)));
                start = None;
            }
            _ => {}
        }
    }
    if let (Some(s), Some(&(m, _))) = (start, pulls.last()) {
        ranges.push((s, m));
    }
    ranges
}

#[cfg(test)]
mod tests {
    use super::*;

    const MSSM: &str = include_str!("../../models/mssm.toml");
    const SPLIT: &str = include_str!("../../models/split_susy.toml");

    #[test]
    fn test_log_grid() {
        let g = log_grid(100.0, 1e6, 5);
        assert_eq!(g.len(), 5);
        assert!((g[0] - 100.0).abs() < 1e-9 && (g[4] - 1e6).abs() < 1e-6);
        assert!((g[1] - 1000.0).abs() < 1e-6);
    }

    #[test]
    fn test_scaling_keeps_hierarchy() {
        let a = Amendment::from_toml(SPLIT).unwrap().with_lightest_mass(5000.0);
        let min = a.added.iter().map(|f| f.mass).fold(f64::MAX, f64::min);
        let max = a.added.iter().map(|f| f.mass).fold(f64::MIN, f64::max);
        assert!((min - 5000.0).abs() < 1e-6);
        assert!((max / min - 1e6).abs() < 1e-3);
    }

    #[test]
    fn test_mssm_best_fit_near_tev() {
        // Two-loop MSSM overshoots alpha_s at M_Z and undershoots at very high
        // masses, so the pull crosses zero in between.
        let a = Amendment::from_toml(MSSM).unwrap();
        let points = scan(&a, 2, &log_grid(100.0, 1e6, 17));
        let (m, pull) = best_fit(&points).unwrap();
        assert_eq!(pull, 0.0);
        assert!(m > 300.0 && m < 1e5, "best fit at {:e} GeV", m);
        let one_sigma = allowed_ranges(&points, 1.0);
        assert_eq!(one_sigma.len(), 1);
        assert!(one_sigma[0].0 < m && m < one_sigma[0].1);
    }
}
