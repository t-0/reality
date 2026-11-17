//! Proton lifetime from gauge-boson exchange at the unification scale.
//!
//! The superheavy X and Y gauge bosons of SU(5)/SO(10) turn quarks into leptons,
//! giving p -> e+ pi0 through a dimension-6 operator with rate
//!   Gamma ~ alpha_GUT^2 m_p^5 / M_X^4.
//! This is the naive dimensional estimate. Hadronic matrix elements, the
//! difference between M_X and the crossing scale, and threshold corrections
//! make it uncertain by about a factor of 10 either way.
//! Supersymmetric models also have dimension-5 decays (p -> K+ nu-bar) that
//! depend on the superpartner spectrum and are not estimated here.

use crate::couplings::constants::M_PROTON;

/// hbar in GeV * s.
const HBAR_GEV_S: f64 = 6.582_119_569e-25;

/// Seconds per year.
const SECONDS_PER_YEAR: f64 = 3.155_76e7;

/// Factor by which the naive estimate is uncertain in each direction.
pub const UNCERTAINTY_FACTOR: f64 = 10.0;

/// Super-Kamiokande lower limit on tau(p -> e+ pi0), 90% CL, in years
/// (Phys. Rev. D 102, 112011, 2020).
pub const SUPER_K_LIMIT_YEARS: f64 = 2.4e34;

/// Approximate Hyper-Kamiokande reach for p -> e+ pi0 after about 10 years
/// of running (90% CL), in years. Hyper-K is due to start around 2027-28.
pub const HYPER_K_REACH_YEARS: f64 = 1.0e35;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    /// Even the longest lifetime in the uncertainty band is below the Super-K limit.
    Excluded,
    /// The uncertainty band straddles the Super-K limit.
    Borderline,
    /// Allowed today, and Hyper-K can see it or rule it out.
    TestableAtHyperK,
    /// Allowed, and beyond Hyper-K's reach.
    OutOfReach,
}

#[derive(Debug, Clone, Copy)]
pub struct ProtonDecay {
    /// Mass of the X/Y gauge bosons, taken as the unification scale (GeV).
    pub m_x: f64,
    pub alpha_gut: f64,
    /// Central estimate of tau(p -> e+ pi0) in years.
    pub lifetime_years: f64,
    pub status: Status,
}

impl ProtonDecay {
    pub fn estimate(m_gut: f64, inverse_alpha_gut: f64) -> Self {
        let alpha_gut = 1.0 / inverse_alpha_gut;
        let width = alpha_gut * alpha_gut * M_PROTON.powi(5) / m_gut.powi(4);
        let lifetime_years = HBAR_GEV_S / width / SECONDS_PER_YEAR;

        let status = classify(lifetime_years, SUPER_K_LIMIT_YEARS, HYPER_K_REACH_YEARS);
        Self {
            m_x: m_gut,
            alpha_gut,
            lifetime_years,
            status,
        }
    }

    pub fn range_years(&self) -> (f64, f64) {
        (
            self.lifetime_years / UNCERTAINTY_FACTOR,
            self.lifetime_years * UNCERTAINTY_FACTOR,
        )
    }
}

/// Place a lifetime, uncertain by UNCERTAINTY_FACTOR either way, against a
/// current limit and a future reach.
pub fn classify(lifetime_years: f64, limit: f64, reach: f64) -> Status {
    let low = lifetime_years / UNCERTAINTY_FACTOR;
    let high = lifetime_years * UNCERTAINTY_FACTOR;
    if high < limit {
        Status::Excluded
    } else if low < limit {
        Status::Borderline
    } else if low < reach {
        Status::TestableAtHyperK
    } else {
        Status::OutOfReach
    }
}

/// Super-Kamiokande lower limit on tau(p -> K+ nu-bar), 90% CL, in years
/// (260 kt yr, arXiv:1408.1195).
pub const SUPER_K_LIMIT_K_NU_YEARS: f64 = 5.9e33;

/// Approximate Hyper-Kamiokande reach for p -> K+ nu-bar, in years.
pub const HYPER_K_REACH_K_NU_YEARS: f64 = 3.0e34;

/// Dimension-5 proton decay p -> K+ nu-bar in supersymmetric GUTs, through
/// colour-triplet higgsino exchange dressed by a higgsino loop.
///
/// Hisano, Kuwahara & Nagata, arXiv:1304.3651, eq. (19), for M_S = mu_H:
///   tau ~ 4e35 yr sin^4(2 beta) (0.1 / A_R)^2 (M_S / 100 TeV)^2 (M_HC / 1e16 GeV)^2
/// The dressing amplitude scales as mu_H / M_S^2, so for a lighter higgsino
/// M_S^2 becomes M_S^4 / mu_H^2. A_R is kept at its reference value 0.1.
/// Uncertain by about a factor of 10 (A_R, hadronic matrix elements, flavour).
#[derive(Debug, Clone, Copy)]
pub struct Dimension5 {
    /// Colour-triplet Higgs mass (GeV).
    pub m_hc: f64,
    /// Sfermion mass scale (GeV).
    pub m_s: f64,
    /// Higgsino mass (GeV), at most m_s.
    pub mu_h: f64,
    pub tan_beta: f64,
    pub lifetime_years: f64,
    pub status: Status,
}

impl Dimension5 {
    pub fn estimate(m_hc: f64, m_s: f64, mu_h: f64, tan_beta: f64) -> Self {
        let mu_h = mu_h.min(m_s);
        let sin_2beta = 2.0 * tan_beta / (1.0 + tan_beta * tan_beta);
        let effective = m_s * m_s / mu_h;
        let lifetime_years = 4e35 * sin_2beta.powi(4) * (effective / 1e5).powi(2) * (m_hc / 1e16).powi(2);
        Self {
            m_hc,
            m_s,
            mu_h,
            tan_beta,
            lifetime_years,
            status: classify(lifetime_years, SUPER_K_LIMIT_K_NU_YEARS, HYPER_K_REACH_K_NU_YEARS),
        }
    }

    pub fn range_years(&self) -> (f64, f64) {
        (
            self.lifetime_years / UNCERTAINTY_FACTOR,
            self.lifetime_years * UNCERTAINTY_FACTOR,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_textbook_value() {
        // M_X = 1e16 GeV, alpha_GUT = 1/25 gives the standard ~1e35 year estimate.
        let p = ProtonDecay::estimate(1e16, 25.0);
        assert!(p.lifetime_years > 1e35 && p.lifetime_years < 3e35, "{:e}", p.lifetime_years);
    }

    #[test]
    fn test_scales_as_fourth_power() {
        let a = ProtonDecay::estimate(1e15, 25.0);
        let b = ProtonDecay::estimate(1e16, 25.0);
        assert!((b.lifetime_years / a.lifetime_years / 1e4 - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_dimension5_reference_point() {
        // tan(beta) = 1 gives sin(2 beta) = 1: the formula's 4e35 years at
        // M_S = mu_H = 100 TeV and M_HC = 1e16 GeV.
        let d = Dimension5::estimate(1e16, 1e5, 1e5, 1.0);
        assert!((d.lifetime_years / 4e35 - 1.0).abs() < 1e-12);
        // Quadratic in M_HC and in M_S (for mu_H = M_S).
        let d2 = Dimension5::estimate(1e17, 1e6, 1e6, 1.0);
        assert!((d2.lifetime_years / 4e39 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_dimension5_kills_tev_susy() {
        // Minimal SUSY SU(5) with TeV sparticles and tan(beta) = 10 is far
        // below the K+ nu-bar limit, the classic result (Murayama & Pierce 2002).
        let d = Dimension5::estimate(1e16, 1e3, 1e3, 10.0);
        assert_eq!(d.status, Status::Excluded);
        // A lighter higgsino suppresses the higgsino dressing.
        let split = Dimension5::estimate(1e16, 1e9, 1e3, 10.0);
        assert!(split.lifetime_years > d.lifetime_years);
    }

    #[test]
    fn test_low_scale_is_excluded() {
        // The SM alpha_1 = alpha_2 crossing near 1e13 GeV: the fate of minimal SU(5).
        assert_eq!(ProtonDecay::estimate(1e13, 42.0).status, Status::Excluded);
        assert_eq!(ProtonDecay::estimate(1e17, 25.0).status, Status::OutOfReach);
    }
}
