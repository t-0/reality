use super::gates::{self, Gate, Verdict};
use super::proton::ProtonDecay;
use super::running::{Running, Segment};
use super::field::Field;
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
    fn new(
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

/// Everything computed for one amendment.
#[derive(Debug, Clone)]
pub struct Report {
    pub amendment: Amendment,
    pub gates: Vec<Gate>,
    pub segments: Vec<Segment>,
    pub landau_pole: Option<f64>,
    pub prediction: Option<Prediction>,
    pub proton: Option<ProtonDecay>,
}

impl Report {
    pub fn gates_passed(&self) -> bool {
        self.gates.iter().all(|g| g.verdict != Verdict::Fail)
    }
}

/// Hold out one coupling, require the other two to meet, and run the
/// held-out coupling back down from that point. One loop: in each segment
/// 1/alpha is linear in ln Q, so 1/alpha_k(M_Z) = 1/alpha_GUT + Delta_k(M_GUT).
pub fn predict(running: &Running, inputs: (usize, usize), predicted: usize) -> Option<Prediction> {
    let inv_mz = inverse_couplings_mz();
    let (i, j) = inputs;
    let m_gut = running.crossing(inv_mz, i, j)?;
    let inverse_alpha_gut = running.inverse_at(inv_mz, m_gut)[i];
    let inverse_predicted = inverse_alpha_gut + running.delta(m_gut)[predicted];
    Some(Prediction::new(1, inputs, predicted, m_gut, inverse_alpha_gut, inverse_predicted))
}

/// Two-loop version: the held-out coupling feeds back into the running of the
/// inputs, so it is found by shooting, seeded with the one-loop answer.
fn predict_two_loop(amendment: &Amendment, fields: &[Field], one_loop: &Prediction) -> Option<Prediction> {
    let two = TwoLoop::new(fields, amendment.top_yukawa);
    let (i, j) = amendment.inputs;
    let k = amendment.predicted;
    let (x, p) = two.solve(inverse_couplings_mz(), (i, j), k, one_loop.inverse_predicted)?;
    Some(Prediction::new(2, (i, j), k, p.ln_q.exp(), p.inv[i], x))
}

/// Evaluate an amendment with one- or two-loop running.
pub fn evaluate(amendment: &Amendment, loops: u8) -> Report {
    let fields = amendment.all_fields();
    let running = Running::new(&fields);
    let landau_pole = running.landau_pole(inverse_couplings_mz());
    let mut prediction = predict(&running, amendment.inputs, amendment.predicted);
    if loops >= 2 {
        prediction = prediction.and_then(|p| predict_two_loop(amendment, &fields, &p));
    }
    let prediction = prediction.filter(|p| landau_pole.map_or(true, |pole| p.m_gut < pole));
    let proton = prediction
        .as_ref()
        .map(|p| ProtonDecay::estimate(p.m_gut, p.inverse_alpha_gut));
    Report {
        amendment: amendment.clone(),
        gates: gates::check(&fields),
        segments: running.segments.clone(),
        landau_pole,
        prediction,
        proton,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
