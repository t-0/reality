use super::gates::{self, Gate, Verdict};
use super::field::{Kind, Su3};
use super::proton::{Dimension5, ProtonDecay};
use super::running::{Running, Segment};
use super::spec::Amendment;
use super::two_loop::TwoLoop;
use crate::couplings::unification::inverse_couplings_mz;

/// Experimental 1-sigma uncertainties on 1/alpha_i(M_Z).
/// From 1/alpha_em(M_Z) = 127.952 +- 0.009, sin^2(theta_W) = 0.23122 +- 0.00004
/// and alpha_s(M_Z) = 0.1180 +- 0.0009 (PDG), propagated to the GUT-normalized couplings.
pub const SIGMA_INV_MZ: [f64; 3] = [0.0052, 0.0055, 0.0646];

/// Relative theory uncertainty on a one-loop prediction of 1/alpha at M_Z.
/// Two-loop running and threshold corrections typically move the predicted
/// alpha_s by 3-8%; 5% is used as a 1-sigma estimate.
pub const ONE_LOOP_THEORY_REL: f64 = 0.05;

/// Relative theory uncertainty at two loops. What remains is mostly unknown
/// threshold corrections at the superpartner and GUT scales, plus three-loop terms.
pub const TWO_LOOP_THEORY_REL: f64 = 0.02;

/// Result of predicting the held-out coupling from the two inputs plus unification.
#[derive(Debug, Clone)]
pub struct Prediction {
    /// Loop order of the running (1 or 2).
    pub loops: u8,
    /// Relative theory uncertainty used in `pull`.
    pub theory_rel: f64,
    pub inputs: (usize, usize),
    pub predicted: usize,
    /// Scale where the two input couplings meet (GeV).
    pub m_gut: f64,
    pub inverse_alpha_gut: f64,
    /// Predicted and measured 1/alpha of the held-out coupling at M_Z.
    pub inverse_predicted: f64,
    pub inverse_measured: f64,
    /// (predicted - measured) / experimental sigma, in 1/alpha.
    pub pull_experimental: f64,
    /// (predicted - measured) / (experimental and one-loop theory sigma in quadrature).
    pub pull: f64,
}

impl Prediction {
    pub(crate) fn new(
        loops: u8,
        inputs: (usize, usize),
        predicted: usize,
        m_gut: f64,
        inverse_alpha_gut: f64,
        inverse_predicted: f64,
    ) -> Self {
        let theory_rel = if loops >= 2 { TWO_LOOP_THEORY_REL } else { ONE_LOOP_THEORY_REL };
        let inverse_measured = inverse_couplings_mz()[predicted];
        let sigma_exp = SIGMA_INV_MZ[predicted];
        let sigma_th = theory_rel * inverse_predicted;
        let diff = inverse_predicted - inverse_measured;
        Self {
            loops,
            theory_rel,
            inputs,
            predicted,
            m_gut,
            inverse_alpha_gut,
            inverse_predicted,
            inverse_measured,
            pull_experimental: diff / sigma_exp,
            pull: diff / (sigma_exp * sigma_exp + sigma_th * sigma_th).sqrt(),
        }
    }

    pub fn alpha_predicted(&self) -> f64 {
        1.0 / self.inverse_predicted
    }
    pub fn alpha_measured(&self) -> f64 {
        1.0 / self.inverse_measured
    }
    pub fn relative_error(&self) -> f64 {
        self.alpha_predicted() / self.alpha_measured() - 1.0
    }
}

/// Spread of the prediction over the corners of the heavy-field ratio ranges.
#[derive(Debug, Clone, Copy)]
pub struct ThresholdBand {
    pub alpha_min: f64,
    pub alpha_max: f64,
    pub m_gut_min: f64,
    pub m_gut_max: f64,
    /// Number of ratio combinations evaluated.
    pub corners: usize,
    /// Range of the dimension-5 lifetime (years) over the corners, for SUSY models.
    pub dimension5_years: Option<(f64, f64)>,
}

impl ThresholdBand {
    pub fn contains(&self, alpha: f64) -> bool {
        self.alpha_min <= alpha && alpha <= self.alpha_max
    }
}

/// Everything computed for one amendment.
#[derive(Debug, Clone)]
pub struct Report {
    pub amendment: Amendment,
    pub gates: Vec<Gate>,
    pub segments: Vec<Segment>,
    pub landau_pole: Option<f64>,
    pub prediction: Option<Prediction>,
    /// Range of the prediction from GUT-scale thresholds, if any are varied.
    pub thresholds: Option<ThresholdBand>,
    pub proton: Option<ProtonDecay>,
    /// p -> K+ nu-bar through colour-triplet higgsinos, for supersymmetric models.
    pub dimension5: Option<Dimension5>,
}

impl Report {
    pub fn gates_passed(&self) -> bool {
        self.gates.iter().all(|g| g.verdict != Verdict::Fail)
    }
}

/// Hold out one coupling, require the other two to meet, and run the
/// held-out coupling back down from that point. One loop: in each segment
/// 1/alpha is linear in ln Q, so 1/alpha_k(M_Z) = 1/alpha_GUT + Delta_k(M_GUT).
///
/// `shift` holds the GUT threshold corrections: unification means
/// 1/alpha_i(M_X) + shift_i = 1/alpha_GUT for all i, with M_X the X boson mass.
pub fn predict(
    running: &Running,
    inputs: (usize, usize),
    predicted: usize,
    shift: [f64; 3],
) -> Option<Prediction> {
    let inv_mz = inverse_couplings_mz();
    let (i, j) = inputs;
    let m_gut = running.crossing(inv_mz, i, j, shift)?;
    let inverse_alpha_gut = running.inverse_at(inv_mz, m_gut)[i] + shift[i];
    let inverse_predicted = inverse_alpha_gut - shift[predicted] + running.delta(m_gut)[predicted];
    Some(Prediction::new(1, inputs, predicted, m_gut, inverse_alpha_gut, inverse_predicted))
}

/// Two-loop version: the held-out coupling feeds back into the running of the
/// inputs, so it is found by shooting, seeded with the one-loop answer.
fn predict_two_loop(
    amendment: &Amendment,
    two: &TwoLoop,
    one_loop: &Prediction,
    shift: [f64; 3],
) -> Option<Prediction> {
    let (i, j) = amendment.inputs;
    let k = amendment.predicted;
    let (x, p) = two.solve(inverse_couplings_mz(), (i, j), k, one_loop.inverse_predicted, shift)?;
    Some(Prediction::new(2, (i, j), k, p.ln_q.exp(), p.inv[i] + shift[i], x))
}

/// Prediction at the given loop order for one set of heavy-field ratios.
fn predict_with(
    amendment: &Amendment,
    running: &Running,
    two: Option<&TwoLoop>,
    ratios: &[f64],
) -> Option<Prediction> {
    let shift = amendment.threshold_shift(ratios);
    let one = predict(running, amendment.inputs, amendment.predicted, shift)?;
    match two {
        Some(t) => predict_two_loop(amendment, t, &one, shift),
        None => Some(one),
    }
}

/// tan(beta) used for dimension-5 proton decay when a model does not set it.
pub const DEFAULT_TAN_BETA: f64 = 10.0;

/// Inputs to the dimension-5 estimate, if the model is supersymmetric at high
/// energy: (sfermion mass M_S, higgsino mass mu_H, index of the colour-triplet
/// Higgs among the heavy fields if declared).
fn dimension5_inputs(amendment: &Amendment, fields: &[super::field::Field]) -> Option<(f64, f64, Option<usize>)> {
    if !super::two_loop::is_supersymmetric(&fields.iter().collect::<Vec<_>>()) {
        return None;
    }
    let squarks = amendment
        .added
        .iter()
        .filter(|f| f.kind == Kind::ComplexScalar && f.su3 != Su3::Singlet)
        .map(|f| f.mass)
        .reduce(f64::max);
    let m_s = squarks.or_else(|| amendment.added.iter().map(|f| f.mass).reduce(f64::max))?;
    let mu_h = amendment
        .added
        .iter()
        .filter(|f| f.kind == Kind::Dirac && f.su3 == Su3::Singlet && f.su2 == 2 && (f.y.abs() - 0.5).abs() < 1e-9)
        .map(|f| f.mass)
        .reduce(f64::min)
        .unwrap_or(m_s);
    let h_c = amendment.heavy.iter().position(|h| {
        let f = &h.field;
        f.kind == Kind::Chiral && matches!(f.su3, Su3::Triplet | Su3::AntiTriplet) && f.su2 == 1 && (f.y.abs() - 1.0 / 3.0).abs() < 1e-9
    });
    Some((m_s, mu_h, h_c))
}

/// Evaluate an amendment with one- or two-loop running.
pub fn evaluate(amendment: &Amendment, loops: u8) -> Report {
    let fields = amendment.all_fields();
    let running = Running::new(&fields);
    let landau_pole = running.landau_pole(inverse_couplings_mz());
    let two = (loops >= 2).then(|| TwoLoop::new(&fields, amendment.top_yukawa));
    let below_pole = |p: &Prediction| landau_pole.map_or(true, |pole| p.m_gut < pole);
    let prediction = predict_with(amendment, &running, two.as_ref(), &amendment.central_ratios())
        .filter(|p| below_pole(p));

    let d5_inputs = dimension5_inputs(amendment, &fields);
    let tan_beta = amendment.tan_beta.unwrap_or(DEFAULT_TAN_BETA);
    let d5 = |p: &Prediction, ratios: &[f64]| {
        d5_inputs.map(|(m_s, mu_h, h_c)| {
            let r = h_c.map_or(1.0, |k| ratios[k]);
            Dimension5::estimate(r * p.m_gut, m_s, mu_h, tan_beta)
        })
    };

    let corner_ratios = amendment.ratio_corners();
    let corners: Vec<(Prediction, Option<f64>)> = corner_ratios
        .iter()
        .filter_map(|r| {
            let p = predict_with(amendment, &running, two.as_ref(), r)?;
            let tau = d5(&p, r).map(|d| d.lifetime_years);
            Some((p, tau))
        })
        .filter(|(p, _)| below_pole(p))
        .collect();
    let thresholds = (!corners.is_empty() && prediction.is_some()).then(|| {
        let fold = |f: fn(&Prediction) -> f64, init: f64, op: fn(f64, f64) -> f64| {
            corners.iter().map(|(p, _)| f(p)).fold(init, op)
        };
        let taus: Vec<f64> = corners.iter().filter_map(|c| c.1).collect();
        ThresholdBand {
            alpha_min: fold(|p| p.alpha_predicted(), f64::MAX, f64::min),
            alpha_max: fold(|p| p.alpha_predicted(), f64::MIN, f64::max),
            m_gut_min: fold(|p| p.m_gut, f64::MAX, f64::min),
            m_gut_max: fold(|p| p.m_gut, f64::MIN, f64::max),
            corners: corners.len(),
            dimension5_years: (!taus.is_empty()).then(|| {
                (taus.iter().cloned().fold(f64::MAX, f64::min), taus.iter().cloned().fold(f64::MIN, f64::max))
            }),
        }
    });
    let proton = prediction
        .as_ref()
        .map(|p| ProtonDecay::estimate(p.m_gut, p.inverse_alpha_gut));
    let dimension5 = prediction.as_ref().and_then(|p| d5(p, &amendment.central_ratios()));
    Report {
        amendment: amendment.clone(),
        gates: gates::check(&fields),
        segments: running.segments.clone(),
        landau_pole,
        prediction,
        thresholds,
        proton,
        dimension5,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::couplings::constants::M_Z;

    fn eval(toml: &str) -> Report {
        evaluate(&Amendment::from_toml(toml).unwrap(), 1)
    }

    #[test]
    fn test_sm_predicts_alpha_s_badly() {
        let p = eval("name = \"SM\"").prediction.unwrap();
        assert!((p.alpha_predicted() - 0.071).abs() < 0.002, "{}", p.alpha_predicted());
        assert!(p.pull_experimental.abs() > 50.0);
        assert!(p.pull.abs() > 5.0, "SM still excluded with theory error: {}", p.pull);
    }

    #[test]
    fn test_inputs_are_reproduced() {
        // Predicting alpha_1 from alpha_2 and alpha_3 in the SM is self-consistent:
        // running alpha_2 to the alpha_2 = alpha_3 crossing and back gives alpha_2.
        let r = eval("name = \"SM\"\ninputs = [\"alpha_2\", \"alpha_3\"]");
        let p = r.prediction.unwrap();
        assert_eq!(p.predicted, 0);
        assert!(p.m_gut > 1e16 && p.m_gut < 1e18);
    }

    #[test]
    fn test_su5_thresholds() {
        let base = Amendment::from_toml(include_str!("../../models/mssm.toml")).unwrap();
        let su5 = Amendment::from_toml(include_str!("../../models/mssm_su5.toml")).unwrap();
        let r0 = evaluate(&base, 2);
        let r = evaluate(&su5, 2);
        // Central ratios are all 1, so the central prediction is unchanged.
        let (a0, a) = (r0.prediction.unwrap().alpha_predicted(), r.prediction.unwrap().alpha_predicted());
        assert!((a - a0).abs() < 1e-6, "{} vs {}", a, a0);
        let band = r.thresholds.unwrap();
        assert_eq!(band.corners, 4);
        assert!(band.alpha_min < a && a < band.alpha_max);

        // A lighter colour-triplet Higgs lowers the predicted alpha_s
        // (the well-known tension with dimension-5 proton decay).
        let shift = |r_hc: f64| su5.threshold_shift(&[1.0, r_hc, 1.0, 1.0]);
        let running = Running::new(&su5.all_fields());
        let light = predict(&running, su5.inputs, su5.predicted, shift(0.1)).unwrap();
        let heavy = predict(&running, su5.inputs, su5.predicted, shift(10.0)).unwrap();
        assert!(light.alpha_predicted() < heavy.alpha_predicted());

        // Hisano-Murayama-Yanagida: at one loop, with all superpartners at M_S,
        // (3/a2 - 2/a3 - 1/a1)(M_Z) = [12/5 ln(M_HC/M_Z) - 2 ln(M_S/M_Z)] / (2 pi),
        // whatever the X and Sigma masses are.
        let inv = inverse_couplings_mz();
        let m_s: f64 = 1000.0;
        for ratios in [[1.0, 0.1, 1.0, 1.0], [1.0, 3.0, 10.0, 10.0], [1.0, 1.0, 0.1, 0.1]] {
            let p = predict(&running, su5.inputs, su5.predicted, su5.threshold_shift(&ratios)).unwrap();
            let lhs = 3.0 * inv[1] - 2.0 * p.inverse_predicted - inv[0];
            let m_hc = ratios[1] * p.m_gut;
            let rhs = (2.4 * (m_hc / M_Z).ln() - 2.0 * (m_s / M_Z).ln()) / (2.0 * std::f64::consts::PI);
            assert!((lhs - rhs).abs() < 1e-6, "ratios {:?}: {} vs {}", ratios, lhs, rhs);
        }
    }

    #[test]
    fn test_no_crossing_gives_no_prediction() {
        // A huge pile of coloured scalars stops alpha_3 from falling, and
        // alpha_3 never meets alpha_2 below Q_MAX.
        let r = eval(
            r#"
            name = "x"
            inputs = ["alpha_2", "alpha_3"]
            [[field]]
            name = "s"
            kind = "complex_scalar"
            su3 = "8"
            copies = 6
            "#,
        );
        assert!(r.prediction.is_none());
    }
}
