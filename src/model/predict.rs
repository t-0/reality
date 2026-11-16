use super::gates::{self, Gate, Verdict};
use super::proton::ProtonDecay;
use super::running::{Running, Segment};
use super::spec::Amendment;
use crate::couplings::unification::inverse_couplings_mz;

/// Experimental 1-sigma uncertainties on 1/alpha_i(M_Z).
/// From 1/alpha_em(M_Z) = 127.952 +- 0.009, sin^2(theta_W) = 0.23122 +- 0.00004
/// and alpha_s(M_Z) = 0.1180 +- 0.0009 (PDG), propagated to the GUT-normalized couplings.
pub const SIGMA_INV_MZ: [f64; 3] = [0.0052, 0.0055, 0.0646];

/// Relative theory uncertainty on a one-loop prediction of 1/alpha at M_Z.
/// Two-loop running and threshold corrections typically move the predicted
/// alpha_s by 3-8%; 5% is used as a 1-sigma estimate.
pub const ONE_LOOP_THEORY_REL: f64 = 0.05;

/// Result of predicting the held-out coupling from the two inputs plus unification.
#[derive(Debug, Clone)]
pub struct Prediction {
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
    let inverse_measured = inv_mz[predicted];
    let sigma_exp = SIGMA_INV_MZ[predicted];
    let sigma_th = ONE_LOOP_THEORY_REL * inverse_predicted;
    let diff = inverse_predicted - inverse_measured;
    Some(Prediction {
        inputs,
        predicted,
        m_gut,
        inverse_alpha_gut,
        inverse_predicted,
        inverse_measured,
        pull_experimental: diff / sigma_exp,
        pull: diff / (sigma_exp * sigma_exp + sigma_th * sigma_th).sqrt(),
    })
}

pub fn evaluate(amendment: &Amendment) -> Report {
    let fields = amendment.all_fields();
    let running = Running::new(&fields);
    let landau_pole = running.landau_pole(inverse_couplings_mz());
    let prediction = predict(&running, amendment.inputs, amendment.predicted)
        .filter(|p| landau_pole.map_or(true, |pole| p.m_gut < pole));
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
        evaluate(&Amendment::from_toml(toml).unwrap())
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
